library(mcmc)
library(DEoptim)
library(sjstats)

#### Estimation functions for the nuclear an cytosolic half life parameters ####

newtotalratio_nucleus <- function(taunu, tp) {
  
  predict_nuc <- 1 - exp(- taunu * tp) 
  
  return(predict_nuc)
}

newtotalratio_cytosol <- function(lambda, taunu, tp) {

  predict_cyt <- 1 - (
   (lambda * exp(-taunu * tp) - taunu * exp(-lambda * tp)) /
     (lambda - taunu)
  )

  return(predict_cyt)
}

negloss <- function(params,
                    tp = timepoints,
                    traf_nuc = transformed_nuc,
                    tot_nuc = total_nuc,
                    traf_cyt = transformed_cyt,
                    tot_cyt = total_cyt) {
  
  taunu <- params[1]
  lambda <- params[2]
  
  
  # reject impossible parameters
  if (any(is.na(c(taunu, lambda)))) {
    return(Inf)
  }
  
  if (any(c(taunu, lambda) <= 0)) {
    return(Inf)
  }
  
  # prediction of the nuclear new/total ratio for the given timepoints
  predict_nuc <- newtotalratio_nucleus(
    taunu = taunu,
    tp = tp
  )
  
  # prediciton of the cytosolic new/total ratio for the given timepoints
  predict_cyt <- newtotalratio_cytosol(
    lambda = lambda,
    taunu = taunu,
    tp = tp
  )

  if (any(is.na(c(predict_nuc, predict_cyt)))) {
    return(Inf)
  } # reject impossible outcomes
  if (any(c(predict_nuc, predict_cyt) < 0)) {
    return(Inf)
  }
  if (any(c(predict_nuc, predict_cyt) > 1)) {
    return(Inf)
  }
  
  # - residual sum of squares for the comparison with the asin(sqrt())
  # transformed, labeling bias corrected, observed new/total ratios
  # in nucleus resp. cytosol
  
  loss_nuc <- sum((asin(sqrt(predict_nuc)) - traf_nuc)^2 * 2 * tot_nuc)
  loss_cyt <- sum((asin(sqrt(predict_cyt)) - traf_cyt)^2 * 2 * tot_cyt)

  res <- - loss_nuc - loss_cyt
  
  return(-res)
}


negloss_onecompartment <- function(param,
                    tp = timepoints,
                    traf_total = transformed,
                    tot_total = total) {
  
  delta <- param
  
  
  # reject impossible parameters
  if (is.na(delta)) {
    return(Inf)
  }
  
  if (delta <= 0) {
    return(Inf)
  }
  
  # prediction of the nuclear new/total ratio for the given timepoints
  predict_total <- newtotalratio_nucleus(
    taunu = delta,
    tp = tp
  )
  
  if (any(is.na(predict_total))) {
    return(Inf)
  } # reject impossible outcomes
  if (any(predict_total < 0)) {
    return(Inf)
  }
  if (any(predict_total > 1)) {
    return(Inf)
  }
  
  # - residual sum of squares for the comparison with the asin(sqrt())
  # transformed, labeling bias corrected, observed new/total ratios
  # in nucleus resp. cytosol
  
  loss_total <- sum((asin(sqrt(predict_total)) - traf_total)^2 * 2 * tot_total)
  
  res <- - loss_total
  
  return(-res)
}

# Implement negloss für nu, using x,y,z parameters

# MCMC Sampling -----------------------------------------------------------------

doMCMC <- function(timepoints,
                   transformed_nuc,
                   total_nuc,
                   transformed_cyt,
                   total_cyt,
                   distinct_samples = 10^5,
                   min_acceptance = 0.15,
                   max_acceptance = 0.45,
                   start_scale = c(1, 1, 1),
                   start_params = c(0.008, 0.008)) {
  # define the loss function in a way that is suitable for the metrop function
  # params is a 3-dim parameter vector containing (tau, r, and lambda)
  loss <- function(params,
                   tp = timepoints,
                   traf_nuc = transformed_nuc,
                   tot_nuc = total_nuc,
                   traf_cyt = transformed_cyt,
                   tot_cyt = total_cyt) {
    res <- negloss(params, tp, traf_nuc, tot_nuc, traf_cyt, tot_cyt)
    
    return(-1 * res)
  }
  # find a suitable scale for the MCMC proposal function
  scale <- start_scale
  res <- metrop(
    loss,
    initial = start_params,
    nbatch = 10^4,
    scale = scale
  )
  
  # increase the scale by a factor of 2 if the acceptance rate is too high
  repeat {
    if (res$accept < max_acceptance) {
      break
    }
    scale <- scale * 2
    res <- metrop(res, nbatch = 10^4, scale = scale)
  }
  # decrease the scale by a factor of 2 if the acceptance rate is too low
  repeat {
    if (res$accept > min_acceptance) {
      break
    }
    scale <- scale / 2
    res <- metrop(res, nbatch = 10^4, scale = scale)
  }
  
  #print("Metrop chain works")
  # perform the final MCMC chain which,
  # in expectation, contains the desired number of distinct samples
  res <- metrop(res,
                nbatch = round(distinct_samples / res$accept),
                scale = scale
  )
  # return the final result
  return(res)
}

doMCMC_onecompartment <- function(timepoints,
                   transformed_total,
                   total_total,
                   distinct_samples = 10^5,
                   min_acceptance = 0.15,
                   max_acceptance = 0.45,
                   start_scale = 1,
                   start_params = 0.008) {
  # define the loss function in a way that is suitable for the metrop function
  # params is a 3-dim parameter vector containing (tau, r, and lambda)
  loss <- function(params,
                   tp = timepoints,
                   traf_total = transformed_total,
                   tot_total = total_total) {
    res <- negloss_onecompartment(params, tp, traf_total, tot_total)
    
    return(-1 * res)
  }
  # find a suitable scale for the MCMC proposal function
  scale <- start_scale
  res <- metrop(
    loss,
    initial = start_params,
    nbatch = 10^4,
    scale = scale
  )
  
  # increase the scale by a factor of 2 if the acceptance rate is too high
  repeat {
    if (res$accept < max_acceptance) {
      break
    }
    scale <- scale * 2
    res <- metrop(res, nbatch = 10^4, scale = scale)
  }
  # decrease the scale by a factor of 2 if the acceptance rate is too low
  repeat {
    if (res$accept > min_acceptance) {
      break
    }
    scale <- scale / 2
    res <- metrop(res, nbatch = 10^4, scale = scale)
  }
  
  #print("Metrop chain works")
  # perform the final MCMC chain which,
  # in expectation, contains the desired number of distinct samples
  res <- metrop(res,
                nbatch = round(distinct_samples / res$accept),
                scale = scale
  )
  # return the final result
  return(res)
}

# Estimation functions ---------------------------------------------------------

# fitTimeseries takes 4 time series and performs an optimization for
# the maximum likelihood parameter estimates,
# first for the nuclear and subsequently for the cytosolic degradation rate
fitTimeseries <- function(timepoints,
                          transformed_nuc,
                          total_nuc,
                          transformed_cyt,
                          total_cyt,
                          include_samples = F,
                          # if TRUE, the whole MCMC chain will also be returned
                          quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                          # quantiles to be reported
                          distinct_samples = 10^4,
                          min_acceptance = 0.15,
                          max_acceptance = 0.45,
                          start_scale = c(1, 1),
                          start_params = c( 0.01, 0.01),
                          genename) {
  # Estimate the maximum likelihood distributions with MCMC sampling and
  # degradation rates with optim
  res <- doMCMC(
    timepoints = timepoints,
    transformed_nuc = transformed_nuc,
    total_nuc = total_nuc,
    transformed_cyt = transformed_cyt,
    total_cyt = total_cyt,
    distinct_samples = distinct_samples,
    min_acceptance = min_acceptance,
    max_acceptance = max_acceptance,
    start_scale = start_scale,
    start_param = c(start_params)
  )
  
  #res$batch <- exp(res$batch)
  
  ML_taunu <- median(res$batch[,1])
  ML_lambda <- median(res$batch[,2])
  #print(c(ML_x, ML_y, ML_z))
  
  est <- tryCatch(
    expr = {
      est <- optim(
        c(ML_taunu, ML_lambda),
        negloss,
        tp = timepoints,
        traf_nuc = transformed_nuc,
        tot_nuc = total_nuc,
        traf_cyt = transformed_cyt,
        tot_cyt = total_cyt,
        method = "Nelder-Mead"
      )$par
    },
    error=function(e){
      cat(genename, "\n")
      est <- optim(
        start_params,
        negloss,
        tp = timepoints,
        traf_nuc = transformed_nuc,
        tot_nuc = total_nuc,
        traf_cyt = transformed_cyt,
        tot_cyt = total_cyt,
        method = "Nelder-Mead"
      )$par
    }
  )
  
  
  taunu <-  est[1]
  lambda <- est[2]
  
  est <- c(taunu, lambda)
  names(est) <- c("taunu", "lambda")
  
  est[c("taunu", "lambda")][which(est[c("taunu", "lambda")] < 0)] <- 0
  
  taunu_samples <- res$batch[, 1]
  lambda_samples <- res$batch[, 2]
  
  # Calculate degradation rate quantiles
  
  quantiles_taunu <- quantile(taunu_samples, probs = quantiles)
  quantile_names <- names(quantiles_taunu)
  names(quantiles_taunu) <- paste("taunu_q", quantile_names, sep = "")
  
  quantiles_lambda <- quantile(lambda_samples, probs = quantiles)
  quantile_names <- names(quantiles_lambda)
  names(quantiles_lambda) <- paste("lambda_q", quantile_names, sep = "")
  
  # Calculate the half lifes and their quantiles (order needs to be reversed!)
  half_life <- log(2) / est[c("taunu", "lambda")]
  names(half_life) <- c(
    "half_life_nuc",
    "half_life_cyt"
  )
  
  quantiles_halflife_nuc <- rev(log(2) / quantiles_taunu)
  names(quantiles_halflife_nuc) <- paste(
    "halflife_nuc_q",
    quantile_names,
    sep = ""
  )
  
  quantiles_halflife_cyt <- rev(log(2) / quantiles_lambda)
  names(quantiles_halflife_cyt) <- paste(
    "halflife_cyt_q",
    quantile_names,
    sep = ""
  )

  
  # Calculate the nuc and cyt Rsquared value (= % of explained variance)
  predict_nuc <- newtotalratio_nucleus(
    taunu = est["taunu"],
    tp = timepoints
  )
  
  rsquared_nuc <- (1 - sum((asin(sqrt(predict_nuc)) - transformed_nuc)^2) / var(transformed_nuc))
  names(rsquared_nuc) <- "Rsquared_nuc"
  
  predict_cyt <- newtotalratio_cytosol(
    lambda = est["lambda"],
    taunu = est["taunu"],
    tp = timepoints
  )
  
  rsquared_cyt <- (1 - sum((asin(sqrt(predict_cyt)) - transformed_cyt)^2) / var(transformed_cyt))
  names(rsquared_cyt) <- "Rsquared_cyt"
  
  # Calculate the mean coverage across the total nuc / total cyt time series
  mean_coverage_nuc <- mean(total_nuc)
  names(mean_coverage_nuc) <- "mean_coverage_nuc"
  mean_coverage_cyt <- mean(total_cyt)
  names(mean_coverage_cyt) <- "mean_coverage_cyt"
  
  # Calculate the coefficient of variation of the total nuc/cyt time series
  cv_nuc <- sd(total_nuc) / mean(total_nuc)
  names(cv_nuc) <- "cv_nuc"
  cv_cyt <- sd(total_cyt) / mean(total_cyt)
  names(cv_cyt) <- "cv_cyt"
  
  # Compose the list of output variables
  final <- list(
    deg_nuc = est["taunu"],
    deg_cyt = est["lambda"],
    logparams = est[c("logtaunu", "loglambda")],
    halflife_nuc = half_life["half_life_nuc"],
    halflife_cyt = half_life["half_life_cyt"],
    quantiles_taunu = quantiles_taunu,
    quantiles_lambda = quantiles_lambda,
    quantiles_halflife_nuc = quantiles_halflife_nuc,
    quantiles_halflife_cyt = quantiles_halflife_cyt,
    mean_coverage_nuc = mean_coverage_nuc,
    mean_coverage_cyt = mean_coverage_cyt,
    cv_nuc = cv_nuc,
    cv_cyt = cv_cyt,
    rsquared_nuc = rsquared_nuc,
    rsquared_cyt = rsquared_cyt,
    nsamples_nuc = nrow(res$batch),
    nsamples_cyt = nrow(res$batch),
    acceptance_nuc = res$accept,
    acceptance_cyt = res$accept,
    time_nuc = res$time[1],
    time_cyt = res$time[1]
  )
  if (include_samples) {
    final$samples_nuc <- res$batch[, 1]
    final$samples_cyt <- res$batch[, 2]
  }
  
  return(final)
}

# fitGene is a wrapper for fitTimeseries
fitGene <- function(gene,
                    # name of the 3'UTR to be fitted
                    timepoints,
                    # vector of effective timepoints
                    rawdata,
                    # data object containing all processed data of one experiment
                    include_samples = F,
                    # if TRUE, the whole MCMC chain will also be returned
                    quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                    # quantiles to be reported
                    distinct_samples = 10^4,
                    min_acceptance = 0.15,
                    max_acceptance = 0.45,
                    start_scale = c(1, 1)) {
  # fetch the required time series corresponding to <gene> from the raw_data object
  # obtain the data for nuclear totals and new/total ratios
  total_nuc <- as.numeric(rawdata$tot_nu[gene, ]) # "uc002sei.1_utr3_0_0_chr2_68405989_r" looks very good
  ratio_nuc <- as.numeric(rawdata$ratio_nu[gene, ])
  transformed_nuc <- as.numeric(asin(sqrt(ratio_nuc)))
  
  # obtain the data for cytosolic totals and new/total ratios
  total_cyt <- as.numeric(rawdata$tot_cy[gene, ])
  ratio_cyt <- as.numeric(rawdata$ratio_cy[gene, ])
  transformed_cyt <- as.numeric(asin(sqrt(ratio_cyt)))
  
  # Find the maximum likelihood seed for MCMC chain, for nuclear degradation
  seed <- DEoptim(
    negloss,
    lower = c(0, 0),
    upper = c(10, 10),
    control = list(itermax = 500, trace = FALSE),
    tp = timepoints,
    traf_nuc = transformed_nuc,
    tot_nuc = total_nuc,
    traf_cyt = transformed_cyt,
    tot_cyt = total_cyt
  )$optim$bestmem
  
  # call fitTimeseries and return its output
  return(
    fitTimeseries(
      timepoints = timepoints,
      transformed_nuc = transformed_nuc,
      total_nuc = total_nuc,
      transformed_cyt = transformed_cyt,
      total_cyt = total_cyt,
      include_samples = include_samples,
      quantiles = quantiles,
      distinct_samples = distinct_samples,
      min_acceptance = min_acceptance,
      max_acceptance = max_acceptance,
      start_scale = start_scale,
      start_params = c(seed),
      genename = gene
    )
  )
}


fitTimeseries_onecompartment <- function(timepoints,
                          transformed_total,
                          total_total,
                          include_samples = F,
                          # if TRUE, the whole MCMC chain will also be returned
                          quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                          # quantiles to be reported
                          distinct_samples = 10^4,
                          min_acceptance = 0.15,
                          max_acceptance = 0.45,
                          start_scale = 1,
                          start_params = 0.01,
                          genename) {
  # Estimate the maximum likelihood distributions with MCMC sampling and
  # degradation rates with optim
  res <- doMCMC_onecompartment(
    timepoints = timepoints,
    transformed_total = transformed_total,
    total_total = total_total,
    distinct_samples = distinct_samples,
    min_acceptance = min_acceptance,
    max_acceptance = max_acceptance,
    start_scale = c(start_scale),
    start_param = c(start_params)
  )
  
  #res$batch <- exp(res$batch)
  
  ML_delta <- median(res$batch[,1])
  
  est <- tryCatch(
    expr = {
      est <- optim(
        ML_delta,
        negloss_onecompartment,
        tp = timepoints,
        traf_total = transformed_total,
        tot_total = total_total,
        method = "Brent",
        lower = 0,
        upper=1
      )$par
    },
    error=function(e){
      cat(genename, "\n")
      est <- optim(
        start_params,
        negloss_onecompartment,
        tp = timepoints,
        traf_total = transformed_total,
        tot_total = total_total,
        method = "Brent",
        lower = 0,
        upper=1
      )$par
    }
  )
  
  
  delta <-  est[1]
  
  est <- delta
  names(est) <- "delta"
  
  est["delta"][which(est["delta"] < 0)] <- 0
  
  delta_samples <- res$batch[, 1]
  
  # Calculate degradation rate quantiles
  
  quantiles_delta <- quantile(delta_samples, probs = quantiles)
  quantile_names <- names(quantiles_delta)
  names(quantiles_delta) <- paste("delta_q", quantile_names, sep = "")
  
  # Calculate the half lifes and their quantiles (order needs to be reversed!)
  half_life <- log(2) / est["delta"]
  names(half_life) <- "half_life_total"
  
  
  quantiles_halflife_total <- rev(log(2) / quantiles_delta)
  names(quantiles_halflife_total) <- paste(
    "halflife_total_q",
    quantile_names,
    sep = ""
  )
  
  
  
  # Calculate the nuc and cyt Rsquared value (= % of explained variance)
  predict_total <- newtotalratio_nucleus(
    taunu = est["delta"],
    tp = timepoints
  )
  
  rsquared_total <- (1 - sum((asin(sqrt(predict_total)) - transformed_total)^2) / var(transformed_total))
  names(rsquared_total) <- "Rsquared_total"

  
  # Calculate the mean coverage across the total nuc / total cyt time series
  mean_coverage_total <- mean(total_total)
  names(mean_coverage_total) <- "mean_coverage_total"
  
  # Calculate the coefficient of variation of the total nuc/cyt time series
  cv_total <- sd(total_total) / mean(total_total)
  names(cv_total) <- "cv_total"
  
  # Compose the list of output variables
  final <- list(
    deg_total = est["delta"],
    logparams = est["logdelta"],
    halflife_total = half_life["half_life_total"],
    quantiles_total = quantiles_delta,
    quantiles_halflife_total = quantiles_halflife_total,
    mean_coverage_total = mean_coverage_total,
    cv_total = cv_total,
    rsquared_total = rsquared_total,
    nsamples_total = nrow(res$batch),
    acceptance_total = res$accept,
    time_total = res$time[1]
  )
  if (include_samples) {
    final$samples_total <- res$batch[, 1]
  }
  
  return(final)
}

# fitGene is a wrapper for fitTimeseries
fitGene_onecompartment <- function(gene,
                    # name of the 3'UTR to be fitted
                    timepoints,
                    # vector of effective timepoints
                    rawdata,
                    # data object containing all processed data of one experiment
                    include_samples = F,
                    # if TRUE, the whole MCMC chain will also be returned
                    quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                    # quantiles to be reported
                    distinct_samples = 10^4,
                    min_acceptance = 0.15,
                    max_acceptance = 0.45,
                    start_scale = 1) {
  # fetch the required time series corresponding to <gene> from the raw_data object
  # obtain the data for nuclear totals and new/total ratios
  total_total <- as.numeric(rawdata$tot_total[gene, ]) 
  ratio_total <- as.numeric(rawdata$ratio_total[gene, ])
  transformed_total <- as.numeric(asin(sqrt(ratio_total)))
  
  # Find the maximum likelihood seed for MCMC chain, for nuclear degradation
  seed <- DEoptim(
    negloss_onecompartment,
    lower = 0,
    upper = 10,
    control = list(itermax = 500, trace = FALSE),
    tp = timepoints,
    traf_total = transformed_total,
    tot_total = total_total
  )$optim$bestmem
  
  # call fitTimeseries and return its output
  return(
    fitTimeseries_onecompartment(
      timepoints = timepoints,
      transformed_total = transformed_total,
      total_total = total_total,
      include_samples = include_samples,
      quantiles = quantiles,
      distinct_samples = distinct_samples,
      min_acceptance = min_acceptance,
      max_acceptance = max_acceptance,
      start_scale = start_scale,
      start_params = c(seed),
      genename = gene
    )
  )
}

# Load data ------------------------------------------------------------------------------------------------------------


load_summarytable <- function(summarytable_file) {

  # *******************************************************************************************************************

  # *** PARAMETERS ***

  # summarytable_file:  file storing summary table

  # *** RETURN ***

  # a data frame containing the summary table (rows: genomic regions, column 1: region name, column 2: measurement
  # description, column 3: library size, column 4: total transcript counts, column 5: labeled transcript counts,
  # column 6: average potential conversion positions, column 7: conversion efficiency column 8: newly synthesized ratio)

  # *******************************************************************************************************************

  df_summarytable <- read.table(summarytable_file, sep = "\t")        # reading in summary table as data frame
  colnames(df_summarytable) <- c("name", "des", "lib", "tot", "mod", "cp", "ce", "nr")    # column names
  df_summarytable <- transform(df_summarytable, name = as.character(name), lib = as.numeric(lib),
                               tot = as.numeric(tot), mod = as.numeric(mod), cp = as.numeric(cp),
                               ce = as.numeric(ce), nr = as.numeric(nr))
  # transforming col 1 (gene names) to strings,
  # 3, 4, 5 (read counts) to intergers and col 6, 7, 8 to floats

  return(df_summarytable)       # returning
}

summarytable_to_rawdata <- function(summarytable_file, timepoints) {

  # *******************************************************************************************************************

  # *** PARAMETERS ***

  # summarytable_file:  file storing summary table
  # timepoints:         vector containing one string representation for each time points measured

  # *** RETURN ***

  # a named list containing 6 matrices storing (in the following order): nuclear modified read counts, nuclear total
  # read counts, nuclear estimated new/total ratios, cytosolic modified read counts, cytosolic total read counts,
  # cytosolic estimated new/total ratios; matrices store time points as columns and genomic region names as rows

  # *******************************************************************************************************************

  stable <- load_summarytable(summarytable_file)    # loading in summary table
  n_timepoints <- length(timepoints)                # number of time points measured
  n_measurements <- n_timepoints * 2                # number of measurements (time series for both nucleus and cytosol)
  n_regions <- nrow(stable) / n_measurements        # number of genomic regions
  region_names <- unique(stable[, "name"])          # names of genomic regions
  n_matrix_entries <- n_timepoints * n_regions      # number of output matrices' entries

  # initializing matrices

  mod_nu <- matrix(rep(0, n_matrix_entries), ncol=n_timepoints)
  colnames(mod_nu) <- timepoints
  rownames(mod_nu) <- region_names
  mod_cy <- matrix(rep(0, n_matrix_entries), ncol=n_timepoints)
  colnames(mod_cy) <- timepoints
  rownames(mod_cy) <- region_names
  tot_nu <- matrix(rep(0, n_matrix_entries), ncol=n_timepoints)
  colnames(tot_nu) <- timepoints
  rownames(tot_nu) <- region_names
  tot_cy <- matrix(rep(0, n_matrix_entries), ncol=n_timepoints)
  colnames(tot_cy) <- timepoints
  rownames(tot_cy) <- region_names
  ratio_nu <- matrix(rep(0, n_matrix_entries), ncol=n_timepoints)
  colnames(ratio_nu) <- timepoints
  rownames(ratio_nu) <- region_names
  ratio_cy <- matrix(rep(0, n_matrix_entries), ncol=n_timepoints)
  colnames(ratio_cy) <- timepoints
  rownames(ratio_cy) <- region_names

  # filling matrices

  for (i in 1:n_regions) {  # iterating through number of genomic regions

    nuc_startidx = (i-1) * n_measurements + 1     # start of nuclear measurements' data rows
    nuc_endidx = nuc_startidx + n_timepoints - 1  # end of nuclear measurements' data rows
    cyt_startidx = nuc_startidx + n_timepoints    # start of cytosolic measurements' data rows
    cyt_endidx = cyt_startidx + n_timepoints - 1  # end of cytosolic measurements' data rows

    mod_nu[i, ] <- stable[nuc_startidx:nuc_endidx, "mod"]
    mod_cy[i, ] <- stable[cyt_startidx:cyt_endidx, "mod"]
    tot_nu[i, ] <- stable[nuc_startidx:nuc_endidx, "tot"]
    tot_cy[i, ] <- stable[cyt_startidx:cyt_endidx, "tot"]
    ratio_nu[i, ] <- stable[nuc_startidx:nuc_endidx, "nr"]
    ratio_cy[i, ] <- stable[cyt_startidx:cyt_endidx, "nr"]

  }

  # returning

  return(list(mod_nu=mod_nu, tot_nu=tot_nu, ratio_nu=ratio_nu, mod_cy=mod_cy, tot_cy=tot_cy, ratio_cy=ratio_cy))
}

summarytable_to_rawdata_onecompartment <- function(summarytable_file, timepoints) {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # summarytable_file:  file storing summary table
  # timepoints:         vector containing one string representation for each time points measured
  
  # *** RETURN ***
  
  # a named list containing 6 matrices storing (in the following order): nuclear modified read counts, nuclear total
  # read counts, nuclear estimated new/total ratios, cytosolic modified read counts, cytosolic total read counts,
  # cytosolic estimated new/total ratios; matrices store time points as columns and genomic region names as rows
  
  # *******************************************************************************************************************
  
  stable <- load_summarytable(summarytable_file)    # loading in summary table
  n_timepoints <- length(timepoints)                # number of time points measured
  n_measurements <- n_timepoints 
  n_regions <- nrow(stable) / n_measurements        # number of genomic regions
  region_names <- unique(stable[, "name"])          # names of genomic regions
  n_matrix_entries <- n_timepoints * n_regions      # number of output matrices' entries
  
  # initializing matrices
  
  mod_total <- matrix(rep(0, n_matrix_entries), ncol=n_timepoints)
  colnames(mod_total) <- timepoints
  rownames(mod_total) <- region_names
  tot_total <- matrix(rep(0, n_matrix_entries), ncol=n_timepoints)
  colnames(tot_total) <- timepoints
  rownames(tot_total) <- region_names
  ratio_total <- matrix(rep(0, n_matrix_entries), ncol=n_timepoints)
  colnames(ratio_total) <- timepoints
  rownames(ratio_total) <- region_names
  
  # filling matrices
  
  for (i in 1:n_regions) {  # iterating through number of genomic regions
    
    total_startidx = (i-1) * n_measurements + 1     # start of nuclear measurements' data rows
    total_endidx = total_startidx + n_timepoints - 1  # end of nuclear measurements' data rows
    
    mod_total[i, ] <- stable[total_startidx:total_endidx, "mod"]
    tot_total[i, ] <- stable[total_startidx:total_endidx, "tot"]
    ratio_total[i, ] <- stable[total_startidx:total_endidx, "nr"]

  }
  
  # returning
  
  return(list(mod_total=mod_total, tot_total=tot_total, ratio_total=ratio_total))
}


estimation_data_table <- function(rawdata,
                                  estimation_results,
                                  elvl_regression,
                                  timepoints,
                                  min_cov,
                                  max_elvl_slope_stringent,
                                  max_elvl_slope_lessstringent,
                                  max_estdev,
                                  max_ciq,
                                  min_r2_stringent,
                                  min_r2_lessstringent) {

  # generate parameter estimation table

  if(any(is.na(estimation_results[,1]))){
    estimation_results <- estimation_results[-which(is.na(estimation_results[,1])),]
  }
  common_loci <- rownames(estimation_results)

  rel_mu <- elvl_regression[common_loci, "avg_elvl_nu"] * estimation_results[common_loci, c("taunu")]
  rel_mu_median <- median(rel_mu)
  rel_mu <- rel_mu / rel_mu_median


  parameter_estimation_table <- cbind(
    estimation_results[, "half_life_nuc"],
    estimation_results[, "half_life_cyt"],
    estimation_results[, "taunu"],
    estimation_results[, "lambda"],
    rel_mu
  )

  colnames(parameter_estimation_table) <- c(
    "half_life_nuc",
    "half_life_cyt",
    "deg_nuc",
    "deg_cyt",
    "rel_mu"
  )

  # generate quantiles table

  quantiles_table <- cbind(
    estimation_results[, "halflife_nuc_q2.5%"],
    estimation_results[, "halflife_nuc_q97.5%"],
    estimation_results[, "halflife_cyt_q2.5%"],
    estimation_results[, "halflife_cyt_q97.5%"],
    estimation_results[, "taunu_q2.5%"],
    estimation_results[, "taunu_q97.5%"],
    estimation_results[, "lambda_q2.5%"],
    estimation_results[, "lambda_q97.5%"]
  )

  colnames(quantiles_table) <- c(
    "half_life_nuc_q2.5%",
    "half_life_nuc_q97.5%",
    "half_life_cyt_q2.5%",
    "half_life_cyt_q97.5%",
    "deg_nuc_q2.5%",
    "deg_nuc_q97.5%",
    "deg_cyt_q2.5%",
    "deg_cyt_q97.5%"
  )

  rownames(quantiles_table) <- common_loci

  # generate metadata table

  negloss_est <- sapply(rownames(estimation_results), function(g) {
    negloss(
      params = estimation_results[g, c("taunu", "lambda")],
      tp = timepoints,
      traf_nuc = asin(sqrt(rawdata$ratio_nu[g, ])),
      tot_nuc = rawdata$tot_nu[g, ],
      traf_cyt = asin(sqrt(rawdata$ratio_cy[g, ])),
      tot_cyt = rawdata$tot_cy[g, ]
    )
  })
  names(negloss_est) <- rownames(estimation_results)

  #elvl_ratio <- (elvl_regression[common_loci, "avg_elvl_cy"] / elvl_regression[common_loci, "avg_elvl_nu"]) / 2
  elvl_ratio <- (elvl_regression[common_loci, "avg_elvl_cy"] / elvl_regression[common_loci, "avg_elvl_nu"])

  metadata_table <- cbind(
    unname(negloss_est),
    estimation_results[, "Rsquared_nuc"],
    estimation_results[, "Rsquared_cyt"],
    estimation_results[, "mean_coverage_nuc"],
    estimation_results[, "mean_coverage_cyt"],
    elvl_regression[common_loci, "avg_elvl_nu"],
    elvl_regression[common_loci, "avg_elvl_cy"],
    elvl_ratio,
    elvl_regression[common_loci, "slope_nu_norm"],
    elvl_regression[common_loci, "slope_cy_norm"]
  )

  colnames(metadata_table) <- c(
    "negloss",
    "r2_nuc",
    "r2_cyt",
    "mean_cov_nuc",
    "mean_cov_cyt",
    "mean_elvl_nuc",
    "mean_elvl_cyt",
    "mean_elvl_ratio_cyt_nuc",
    "expr_lvl_slope_nuc",
    "expr_lvl_slope_cyt"
  )

  # generate reliability table

  reliability_booleans <- cbind(
    as.numeric(metadata_table[, "mean_cov_nuc"] >= min_cov),
    as.numeric(metadata_table[, "mean_cov_cyt"] >= min_cov),
    as.numeric(((parameter_estimation_table[, "deg_nuc"] - quantiles_table[, "deg_nuc_q2.5%"]) /
                  parameter_estimation_table[, "deg_nuc"]) <= max_ciq),
    as.numeric(((quantiles_table[, "deg_nuc_q97.5%"] - parameter_estimation_table[, "deg_nuc"]) /
                  parameter_estimation_table[, "deg_nuc"]) <= max_ciq),
    as.numeric(((parameter_estimation_table[, "deg_cyt"] - quantiles_table[, "deg_cyt_q2.5%"]) /
                  parameter_estimation_table[, "deg_cyt"]) <= max_ciq),
    as.numeric(((quantiles_table[, "deg_cyt_q97.5%"] - parameter_estimation_table[, "deg_cyt"]) /
                  parameter_estimation_table[, "deg_cyt"]) <= max_ciq),
    as.numeric(abs(metadata_table[, "expr_lvl_slope_nuc"]) <= max_elvl_slope_stringent),
    as.numeric(abs(metadata_table[, "expr_lvl_slope_cyt"]) <= max_elvl_slope_stringent),
    as.numeric(abs(metadata_table[, "expr_lvl_slope_nuc"]) <= max_elvl_slope_lessstringent),
    as.numeric(abs(metadata_table[, "expr_lvl_slope_cyt"]) <= max_elvl_slope_lessstringent),
    as.numeric(estimation_results[, "Rsquared_nuc"] >= min_r2_stringent),
    as.numeric(estimation_results[, "Rsquared_nuc"] >= min_r2_lessstringent),
    as.numeric(estimation_results[, "Rsquared_cyt"] >= min_r2_stringent),
    as.numeric(estimation_results[, "Rsquared_cyt"] >= min_r2_lessstringent)
  )

  rownames(reliability_booleans) <- rownames(estimation_results)
  colnames(reliability_booleans) <- c(
    "relab_cov_nuc",
    "relab_cov_cyt",
    "relab_ciq_nuc_lower",
    "relab_ciq_nuc_upper",
    "relab_ciq_cyt_lower",
    "relab_ciq_cyt_upper",
    "relab_elvlslope_nuc",
    "relab_elvlslope_cyt",
    "relab_less_elvlslope_nuc",
    "relab_less_elvlslope_cyt",
    "relab_r2_nuc",
    "relab_less_r2_nuc",
    "relab_r2_cyt",
    "relab_less_r2_cyt"
  )

  reliability_scores <- cbind(
    rowSums(reliability_booleans[
      common_loci,
      c(
        "relab_cov_nuc", "relab_elvlslope_nuc",
        "relab_ciq_nuc_lower", "relab_ciq_nuc_upper",
        "relab_r2_nuc"
      )
    ]),
    rowSums(reliability_booleans[
      common_loci,
      c(
        "relab_cov_nuc", "relab_less_elvlslope_nuc",
        "relab_ciq_nuc_lower", "relab_ciq_nuc_upper",
        "relab_less_r2_nuc"
      )
    ]),
    rowSums(reliability_booleans[
      common_loci,
      c(
        "relab_cov_nuc",
        "relab_cov_cyt",
        "relab_ciq_nuc_lower", "relab_ciq_nuc_upper",
        "relab_ciq_cyt_lower", "relab_ciq_cyt_upper",
        "relab_r2_nuc",
        "relab_r2_cyt",
        "relab_elvlslope_nuc", "relab_elvlslope_cyt"
      )
    ]),
    rowSums(reliability_booleans[
      common_loci,
      c(
        "relab_cov_nuc",
        "relab_cov_cyt",
        "relab_ciq_nuc_lower", "relab_ciq_nuc_upper",
        "relab_ciq_cyt_lower", "relab_ciq_cyt_upper",
        "relab_less_r2_nuc",
        "relab_less_r2_cyt",
        "relab_less_elvlslope_nuc", "relab_less_elvlslope_cyt"
      )
    ])
  )
  rownames(reliability_scores) <- common_loci
  colnames(reliability_scores) <- c(
    "relab_score_stringent_nuc", "relab_score_lessstringent_nuc",
    "relab_score_stringent_both", "relab_score_lessstringent_both"
  )

  reliability_outcomes <- cbind(
    as.numeric(reliability_scores[common_loci, "relab_score_stringent_nuc"] == 5),
    as.numeric(reliability_scores[common_loci, "relab_score_lessstringent_nuc"] == 5),
    as.numeric(reliability_scores[common_loci, "relab_score_stringent_both"] == 10),
    as.numeric(reliability_scores[common_loci, "relab_score_lessstringent_both"] == 10)
  )
  rownames(reliability_outcomes) <- common_loci
  colnames(reliability_outcomes) <- c(
    "reliability_stringent_nuc", "reliability_lessstringent_nuc",
    "reliability_stringent_both", "reliability_lessstringent_both"
  )

  reliability_table <- cbind(
    reliability_booleans[common_loci, ], reliability_scores[common_loci, ],
    reliability_outcomes[common_loci, ]
  )

  # all data table

  all_data_table <- cbind(
    parameter_estimation_table[common_loci, ], quantiles_table[common_loci, ],
    metadata_table[common_loci, ], reliability_table[common_loci, ]
  )

  # returning

  return(list(
    params = parameter_estimation_table, quants = quantiles_table,
    meta = metadata_table, relab = reliability_table, all = all_data_table
  ))
}

estimation_data_table_onecompartment <- function(rawdata,
                                  estimation_results,
                                  elvl_regression,
                                  timepoints,
                                  min_cov,
                                  max_elvl_slope_stringent,
                                  max_elvl_slope_lessstringent,
                                  max_estdev,
                                  max_ciq,
                                  min_r2_stringent,
                                  min_r2_lessstringent) {
  
  # generate parameter estimation table
  
  if(any(is.na(estimation_results[,1]))){
    estimation_results <- estimation_results[-which(is.na(estimation_results[,1])),]
  }
  common_loci <- rownames(estimation_results)
  
  rel_mu <- elvl_regression[common_loci, "avg_elvl_total"] * estimation_results[common_loci, c("delta")]
  rel_mu_median <- median(rel_mu)
  rel_mu <- rel_mu / rel_mu_median
  
  
  parameter_estimation_table <- cbind(
    estimation_results[, "half_life_total"],
    estimation_results[, "delta"],
    rel_mu
  )
  
  colnames(parameter_estimation_table) <- c(
    "half_life_total",
    "deg_total",
    "rel_mu"
  )
  
  # generate quantiles table
  
  quantiles_table <- cbind(
    estimation_results[, "halflife_total_q2.5%"],
    estimation_results[, "halflife_total_q97.5%"],
    estimation_results[, "delta_q2.5%"],
    estimation_results[, "delta_q97.5%"]
  )
  
  colnames(quantiles_table) <- c(
    "half_life_total_q2.5%",
    "half_life_total_q97.5%",
    "deg_total_q2.5%",
    "deg_total_q97.5%"
  )
  
  rownames(quantiles_table) <- common_loci
  
  # generate metadata table
  
  negloss_est <- sapply(rownames(estimation_results), function(g) {
    negloss_onecompartment(
      param = estimation_results[g, "delta"],
      tp = timepoints,
      traf_total = asin(sqrt(rawdata$ratio_total[g, ])),
      tot_total = rawdata$tot_nu[g, ]
    )
  })
  names(negloss_est) <- rownames(estimation_results)
  
  metadata_table <- cbind(
    unname(negloss_est),
    estimation_results[, "Rsquared_total"],
    estimation_results[, "mean_coverage_total"],
    elvl_regression[common_loci, "avg_elvl_total"],
    elvl_regression[common_loci, "slope_total_norm"]
  )
  
  colnames(metadata_table) <- c(
    "negloss",
    "r2_total",
    "mean_cov_total",
    "mean_elvl_total",
    "expr_lvl_slope_total"
  )
  
  # generate reliability table
  
  reliability_booleans <- cbind(
    as.numeric(metadata_table[, "mean_cov_total"] >= min_cov),
    as.numeric(((parameter_estimation_table[, "deg_total"] - quantiles_table[, "deg_total_q2.5%"]) /
                  parameter_estimation_table[, "deg_total"]) <= max_ciq),
    as.numeric(((quantiles_table[, "deg_total_q97.5%"] - parameter_estimation_table[, "deg_total"]) /
                  parameter_estimation_table[, "deg_total"]) <= max_ciq),
    as.numeric(abs(metadata_table[, "expr_lvl_slope_total"]) <= max_elvl_slope_stringent),
    as.numeric(abs(metadata_table[, "expr_lvl_slope_total"]) <= max_elvl_slope_lessstringent),
    as.numeric(estimation_results[, "Rsquared_total"] >= min_r2_stringent),
    as.numeric(estimation_results[, "Rsquared_total"] >= min_r2_lessstringent)
  )
  
  rownames(reliability_booleans) <- rownames(estimation_results)
  colnames(reliability_booleans) <- c(
    "relab_cov_total",
    "relab_ciq_total_lower",
    "relab_ciq_total_upper",
    "relab_elvlslope_total",
    "relab_less_elvlslope_total",
    "relab_r2_total",
    "relab_less_r2_total"
  )
  
  reliability_scores <- cbind(
    rowSums(reliability_booleans[
      common_loci,
      c(
        "relab_cov_total", "relab_elvlslope_total",
        "relab_ciq_total_lower", "relab_ciq_total_upper",
        "relab_r2_total"
      )
    ]),
    rowSums(reliability_booleans[
      common_loci,
      c(
        "relab_cov_total", "relab_less_elvlslope_total",
        "relab_ciq_total_lower", "relab_ciq_total_upper",
        "relab_less_r2_total"
      )
    ])
  )
  rownames(reliability_scores) <- common_loci
  colnames(reliability_scores) <- c(
    "relab_score_stringent_total", "relab_score_lessstringent_total"
  )
  
  reliability_outcomes <- cbind(
    as.numeric(reliability_scores[common_loci, "relab_score_stringent_total"] == 5),
    as.numeric(reliability_scores[common_loci, "relab_score_lessstringent_total"] == 5)
  )
  rownames(reliability_outcomes) <- common_loci
  colnames(reliability_outcomes) <- c(
    "reliability_stringent_total", "reliability_lessstringent_total"
  )
  
  reliability_table <- cbind(
    reliability_booleans[common_loci, ], reliability_scores[common_loci, ],
    reliability_outcomes[common_loci, ]
  )
  
  # all data table
  
  all_data_table <- cbind(
    parameter_estimation_table[common_loci, ], quantiles_table[common_loci, ],
    metadata_table[common_loci, ], reliability_table[common_loci, ]
  )
  
  # returning
  
  return(list(
    params = parameter_estimation_table, quants = quantiles_table,
    meta = metadata_table, relab = reliability_table, all = all_data_table
  ))
}

### Data Processing ###################################################################################################


expr_robust_average <- function(summary_table, n_timepoints) {

  # *** PARAMETERS ***

  # summary_table:      summary table matrix
  # n_timepoints:       number of time points at which measurements were taken

  # *** NOTE ***

  # summary table must be sorted by gene names and time points (nuclear time series, cytosolic time series)

  # *** RETURN ***

  # a vector containing the average of the total counts distribution, computed for the middle 50% of the distribution
  # (for robustness), for each measurement time point

  n_rows <- nrow(summary_table) # number of data rows in summary table
  n_genes <- n_rows / (2 * n_timepoints) # number of genes in the summary table
  quant_25 <- ceiling(n_genes / 100 * 25) # 25% quantile index
  quant_75 <- floor(n_genes / 100 * 75) # 75% quantile index
  interval_50_size <- quant_75 - quant_25 + 1 # 50% interval sample size
  total_collections <- sapply(1:(n_timepoints * 2), function(i) { # collecting total counts for each time point
    summary_table[seq(i, n_rows, by = (n_timepoints * 2)), "tot"]
  })
  total_avgs <- sapply(1:ncol(total_collections), function(i) { # computing averages of time points' total counts
    total_coll <- total_collections[, i]
    total_coll <- sort(total_coll)
    return(sum(total_coll[quant_25:quant_75]) / interval_50_size)
  })
  return(total_avgs)
}

expr_robust_average_onecompartment <- function(summary_table, n_timepoints) {
  
  # *** PARAMETERS ***
  
  # summary_table:      summary table matrix
  # n_timepoints:       number of time points at which measurements were taken
  
  # *** NOTE ***
  
  # summary table must be sorted by gene names and time points
  
  # *** RETURN ***
  
  # a vector containing the average of the total counts distribution, computed for the middle 50% of the distribution
  # (for robustness), for each measurement time point
  
  n_rows <- nrow(summary_table) # number of data rows in summary table
  n_genes <- n_rows / (n_timepoints) # number of genes in the summary table
  quant_25 <- ceiling(n_genes / 100 * 25) # 25% quantile index
  quant_75 <- floor(n_genes / 100 * 75) # 75% quantile index
  interval_50_size <- quant_75 - quant_25 + 1 # 50% interval sample size
  total_collections <- sapply(1:(n_timepoints), function(i) { # collecting total counts for each time point
    summary_table[seq(i, n_rows, by = (n_timepoints)), "tot"]
  })
  total_avgs <- sapply(1:ncol(total_collections), function(i) { # computing averages of time points' total counts
    total_coll <- total_collections[, i]
    total_coll <- sort(total_coll)
    return(sum(total_coll[quant_25:quant_75]) / interval_50_size)
  })
  return(total_avgs)
}

expr_level_regression <- function(summary_table, time_series, norm = "lib") {

  # *******************************************************************************************************************

  # *** PARAMETERS ***

  # summary_table:      summary table file
  # time_series:        vector of time points at which measurements were taken (not including time point 0)
  # norm:               value by which total counts are normalized to compute expression levels; choose 'lib' to use
  #                     the library size, choose 'avg' to use the average of the total counts distribution, computed
  #                     for the middle 50% of the distribution (for robustness) (default: "lib")

  # *** NOTE ***

  # summary table must be sorted by gene names and time points (nuclear time series, cytosolic time series)

  # *** RETURN ***

  # a matrix containing each gene's linear regression model (rows: genes, columns (for both nucleus and cytosol, in
  # that order): intercept, slope, slope normalized w.r.t average expression level, intercept p-value, slope p-value,
  # model p-value, R2 value, coefficient of variation, average expression level)

  # *******************************************************************************************************************

  summary_table <- load_summarytable(summary_table) # loading in summary tables
  n_timepoints <- length(time_series) # number of measurement time points
  gene_names <- unique(summary_table[, "name"]) # gene names

  # computing expression levels (rows: genes, columns: time series' expression levels)
  if (norm == "lib") { # total/library_size ratios
    elvl <- summary_table[, "tot"] / summary_table[, "lib"]
  }
  else if (norm == "avg") { # total/total_avg ratios
    total_avgs <- expr_robust_average(summary_table, n_timepoints)
    elvl <- summary_table[, "tot"] / total_avgs
  }

  elvl <- matrix(elvl, ncol = 2 * n_timepoints, byrow = TRUE) # matrix for expression
  cy_elvl <- elvl[, 1:n_timepoints] # nuclear expression levels
  nu_elvl <- elvl[, (n_timepoints + 1):(2 * n_timepoints)] # cytosolic expression levels

  # iterating through all genes, performing linear regression (matrix; columns: genes, rows (for both nucleus and
  # cytosol, in that order): intercept, slope, intercept p-value, slope p-value, model p-value, R2 value
  lm_matrix <- sapply(1:nrow(elvl), function(i) {

    # data frames storing time series and corresponding expression levels of current gene
    gene_nu_elvl <- data.frame(matrix(c(time_series, nu_elvl[i, ]), ncol = 2))
    colnames(gene_nu_elvl) <- c("time", "elvl")
    gene_cy_elvl <- data.frame(matrix(c(time_series, cy_elvl[i, ]), ncol = 2))
    colnames(gene_cy_elvl) <- c("time", "elvl")

    # computing average expression levels
    nu_elvl_avg <- sum(gene_nu_elvl$elvl) / n_timepoints
    cy_elvl_avg <- sum(gene_cy_elvl$elvl) / n_timepoints

    # fitting linear model
    lm_nu <- lm(elvl ~ time, data = gene_nu_elvl)
    lm_cy <- lm(elvl ~ time, data = gene_cy_elvl)

    # computing model p-values
    lm_summary_nu <- summary(lm_nu) # contains coefficients and corresponding p-values as well as F-statistic
    lm_summary_cy <- summary(lm_cy)
    lm_nu_cv <- cv(lm_nu) # model coefficient of variation
    lm_cy_cv <- cv(lm_cy)

    fstat_nu <- lm_summary_nu$fstatistic # getting F-statistic
    model_p_nu <- pf(fstat_nu[1], fstat_nu[2], fstat_nu[3], lower = FALSE)[[1]] # getting model p-value
    fstat_cy <- lm_summary_cy$fstatistic
    model_p_cy <- pf(fstat_cy[1], fstat_cy[2], fstat_cy[3], lower = FALSE)[[1]]

    # creating return column
    intercept_nu <- lm_summary_nu$coefficients[1, 1]
    slope_nu <- lm_summary_nu$coefficients[2, 1]
    slope_nu_norm <- slope_nu / nu_elvl_avg
    intercept_cy <- lm_summary_cy$coefficients[1, 1]
    slope_cy <- lm_summary_cy$coefficients[2, 1]
    slope_cy_norm <- slope_cy / cy_elvl_avg
    intercept_p_nu <- lm_summary_nu$coefficients[1, 4]
    slope_p_nu <- lm_summary_nu$coefficients[2, 4]
    intercept_p_cy <- lm_summary_cy$coefficients[1, 4]
    slope_p_cy <- lm_summary_cy$coefficients[2, 4]

    return_col <-
      c(
        intercept_nu, slope_nu, slope_nu_norm, intercept_p_nu, slope_p_nu, model_p_nu,
        lm_summary_nu$r.squared, lm_nu_cv, nu_elvl_avg,
        intercept_cy, slope_cy, slope_cy_norm, intercept_p_cy, slope_p_cy, model_p_cy,
        lm_summary_cy$r.squared, lm_cy_cv, cy_elvl_avg
      )
    return(return_col)
  })

  # naming rows and columns of output matrix and returning

  colnames(lm_matrix) <- gene_names
  rownames(lm_matrix) <-
    c(
      "intercept_nu", "slope_nu", "slope_nu_norm", "intercept_p_nu", "slope_p_nu", "model_p_nu",
      "R2_nu", "cv_nu", "avg_elvl_nu",
      "intercept_cy", "slope_cy", "slope_cy_norm", "intercept_p_cy", "slope_p_cy", "model_p_cy",
      "R2_cy", "cv_cy", "avg_elvl_cy"
    )
  lm_matrix <- t(lm_matrix)
  return(lm_matrix)
}

expr_level_regression_onecompartment <- function(summary_table, time_series, norm = "lib") {
  
  # *******************************************************************************************************************
  
  # *** PARAMETERS ***
  
  # summary_table:      summary table file
  # time_series:        vector of time points at which measurements were taken (not including time point 0)
  # norm:               value by which total counts are normalized to compute expression levels; choose 'lib' to use
  #                     the library size, choose 'avg' to use the average of the total counts distribution, computed
  #                     for the middle 50% of the distribution (for robustness) (default: "lib")
  
  # *** NOTE ***
  
  # summary table must be sorted by gene names and time points (nuclear time series, cytosolic time series)
  
  # *** RETURN ***
  
  # a matrix containing each gene's linear regression model (rows: genes, columns (for both nucleus and cytosol, in
  # that order): intercept, slope, slope normalized w.r.t average expression level, intercept p-value, slope p-value,
  # model p-value, R2 value, coefficient of variation, average expression level)
  
  # *******************************************************************************************************************
  
  summary_table <- load_summarytable(summary_table) # loading in summary tables
  n_timepoints <- length(time_series) # number of measurement time points
  gene_names <- unique(summary_table[, "name"]) # gene names
  
  # computing expression levels (rows: genes, columns: time series' expression levels)
  if (norm == "lib") { # total/library_size ratios
    elvl <- summary_table[, "tot"] / summary_table[, "lib"]
  }
  if (norm == "avg") { # total/total_avg ratios
    total_avgs <- expr_robust_average(summary_table, n_timepoints)
    elvl <- summary_table[, "tot"] / total_avgs
  }
  
  elvl <- matrix(elvl, ncol = n_timepoints, byrow = TRUE) # matrix for expression
  total_elvl <- elvl[, 1:n_timepoints] # total expression levels
  
  # iterating through all genes, performing linear regression (matrix; columns: genes, rows (for both nucleus and
  # cytosol, in that order): intercept, slope, intercept p-value, slope p-value, model p-value, R2 value
  lm_matrix <- sapply(1:nrow(elvl), function(i) {
    
    # data frames storing time series and corresponding expression levels of current gene
    gene_total_elvl <- data.frame(matrix(c(time_series, total_elvl[i, ]), ncol = 2))
    colnames(gene_total_elvl) <- c("time", "elvl")
    
    # computing average expression levels
    total_elvl_avg <- sum(gene_total_elvl$elvl) / n_timepoints

    # fitting linear model
    lm_total <- lm(elvl ~ time, data = gene_total_elvl)

    # computing model p-values
    lm_summary_total <- summary(lm_total) # contains coefficients and corresponding p-values as well as F-statistic
    lm_total_cv <- cv(lm_total) # model coefficient of variation

    fstat_total <- lm_summary_total$fstatistic # getting F-statistic
    model_p_total <- pf(fstat_total[1], fstat_total[2], fstat_total[3], lower = FALSE)[[1]] # getting model p-value
    
    
    # creating return column
    intercept_total <- lm_summary_total$coefficients[1, 1]
    slope_total <- lm_summary_total$coefficients[2, 1]
    slope_total_norm <- slope_total / total_elvl_avg

    intercept_p_total <- lm_summary_total$coefficients[1, 4]
    slope_p_total <- lm_summary_total$coefficients[2, 4]
    
    return_col <-
      c(
        intercept_total, slope_total, slope_total_norm, intercept_p_total, slope_p_total, model_p_total,
        lm_summary_total$r.squared, lm_total_cv, total_elvl_avg
      )
    return(return_col)
  })
  
  # naming rows and columns of output matrix and returning
  
  colnames(lm_matrix) <- gene_names
  rownames(lm_matrix) <-
    c(
      "intercept_total", "slope_total", "slope_total_norm", "intercept_p_total", "slope_p_total", "model_p_total",
      "R2_total", "cv_total", "avg_elvl_total"
    )
  lm_matrix <- t(lm_matrix)
  return(lm_matrix)
}
