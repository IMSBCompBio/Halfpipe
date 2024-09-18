library(argparse)
library(doParallel)
library(foreach)

set.seed(1234) # for MCMC sampling reproducibility

parser <- ArgumentParser(description='Process some integers')
parser$add_argument('-i','--input', help='<Required> Set flag', required=TRUE, type="character")
parser$add_argument('-o','--output', help='<Required> Set flag', required=TRUE, type="character")
parser$add_argument('-tp','--timepoints', nargs='+', help='<Required> Set flag', required=TRUE, type="integer")
parser$add_argument('-m','--model', help='<Required> Set flag', required=TRUE, default="onecompartment", type="character")
parser$add_argument('-c','--cores', help='<Required> Set flag', required=FALSE, default=1, type="integer")



args <- parser$parse_args()

# Load functions ---------------------------------------------------------------

source("src/fit_functions.R")

# Load Data --------------------------------------------------------------------

print(args$timepoints)

if (args$model == "onecompartment") {
  
  cat("Loading Raw Data:\n")
  rawdata <- summarytable_to_rawdata_onecompartment(
    summarytable_file = args$input,
    timepoints = args$timepoints
  )
  cat("Data loaded.\n")
  genes <- rownames(rawdata$tot_total)[1:5]
  
  # Estimation -----------------------------------------------------------------
  
  parallel_clusters <- makePSOCKcluster(args$cores, outfile = "")
  registerDoParallel(parallel_clusters)
  
  # fit 3'UTRs in both experiments
  cat("Fitting Parameters:\n")
  
  results <- foreach(
    j = seq_along(genes),
    .packages = c("mcmc", "DEoptim"),
    .export = c("fitGene"),
    .verbose = FALSE
  ) %dopar% { # run estimation on multiple cores
    
    if (j %% 10 == 0) {
      cat(j, " , ")
    }
    gene <- genes[j]
    
    
    fit <- tryCatch(
      {fit <- fitGene_onecompartment(
        gene,
        c(15,30,45),
        rawdata,
        distinct_samples = 10^4,
        include_samples = T
      )},
      error=function(e){
        cat(paste("UTR ", j, "failed", sep="\t"), sep="\n", append = TRUE)
        return(list(NA))
      },
      warning=function(w){
        cat(paste("UTR ", j, "failed", sep="\t"), sep="\n", append = TRUE)
        return(list(NA))
      }
    )
    
    # cumulative function for MCMC samples
    
    if (length(fit) > 1){
      cumf_total<- ecdf(fit$samples_total)
      if (j > 100) {
        fit$samples_total <- NULL
      }
      
      # within-sample quantiles
      within_quantiles <- cumf_total(fit$deg_total)
      
      names(within_quantiles) <- "total_quant"
      # returning
      fit$within_quant <- within_quantiles
    }
    
    
    
    return(fit)
  }
  names(results) <- genes
  cat("\n")
  
  stopCluster(parallel_clusters) # giving free parallelization cores
  
  # Create and store summary tables tables -------------------------------------
  
  cat("Creating Summary Tables...\n")
  results_table <- NULL
  for (j in seq_along(results)) {
    entry <- results[[j]]
    if (length(entry) <= 1){
      entry_1 <- rep(NA, 16)
      entry_2 <- rep(NA, 16)
    } else {
      entry <- with(
        entry,
        c(
          halflife_total,
          deg_total,
          quantiles_halflife_total,
          quantiles_total,
          within_quant,
          mean_coverage_total,
          cv_total,
          rsquared_total
        )
      )
    }
    
    results_table <- rbind(results_table, entry)
  }
  
  column_names <- colnames(results_table)
  results_table <- data.frame(results_table)
  colnames(results_table) <- column_names
  rownames(results_table) <- genes
  print(results_table)
  
  # Create Half-Life table with reliability criteria ---------------------------
  
  elvl_regression <- expr_level_regression_onecompartment(
    args$input,
    time_series = args$timepoints,
    norm = "avg"
  )
  
  half_life_table <- estimation_data_table_onecompartment(
    rawdata = rawdata,
    estimation_results = results_table,
    elvl_regression = elvl_regression,
    timepoints = args$timepoints,
    min_cov = 30,
    max_elvl_slope_stringent = 0.002,
    max_elvl_slope_lessstringent = 0.0025,
    max_estdev = 0.33,
    max_ciq = 0.3,
    min_r2_stringent = 0.5,
    min_r2_lessstringent = 0.4
  )$all
  
  cat("Saving Results...\n")
  # Save the results
  
  save(
    results,
    results_table,
    half_life_table,
    file = args$output
  )
  cat("Computation Done.")
  
}

if (args$model == "twocompartment") {
  
  cat("Loading Raw Data:\n")
  rawdata <- summarytable_to_rawdata(
    summarytable_file = args$input,
    timepoints = args$timepoints
  )
  cat("Data loaded.\n")
  genes <- rownames(rawdata$tot_cy)[1:5]
  
  # Estimation -----------------------------------------------------------------
  
  parallel_clusters <- makePSOCKcluster(args$cores, outfile = "")
  registerDoParallel(parallel_clusters)
  
  # fit 3'UTRs in both experiments
  cat("Fitting Parameters:\n")
  
  results <- foreach(
    j = seq_along(genes),
    .packages = c("mcmc", "DEoptim"),
    .export = c("fitGene"),
    .verbose = FALSE
  ) %dopar% { # run estimation on multiple cores
    
    if (j %% 10 == 0) {
      cat(j, " , ")
    }
    gene <- genes[j]
    
    
    fit <- tryCatch(
      {fit <- fitGene(
        gene,
        args$timepoints,
        rawdata,
        distinct_samples = 10^4,
        include_samples = T
      )},
      error=function(e){
        cat(paste("UTR ", j, "failed", sep="\t"), sep="\n", append = TRUE)
        return(list(NA))
      },
      warning=function(w){
        cat(paste("UTR ", j, "failed", sep="\t"), sep="\n", append = TRUE)
        return(list(NA))
      }
    )
    
    # cumulative function for MCMC samples
    
    if (length(fit) > 1){
      cumf_nuc<- ecdf(fit$samples_nuc)
      cumf_cyt <- ecdf(fit$samples_cyt)
      if (j > 100) {
        fit$samples_nuc <- NULL
        fit$samples_cyt <- NULL
      }
      
      # within-sample quantiles
      within_quantiles <- c(
        cumf_nuc(fit$deg_nuc),
        cumf_cyt(fit$deg_cyt)
      )
      names(within_quantiles) <- c(
        "nuc_quant",
        "cyt_quant"
      )
      # returning
      fit$within_quant <- within_quantiles
    }
    
    
    
    return(fit)
  }
  names(results) <- genes
  cat("\n")
  
  stopCluster(parallel_clusters) # giving free parallelization cores
  
  # Create and store summary tables tables -------------------------------------
  
  cat("Creating Summary Tables...\n")
  results_table <- NULL
  for (j in seq_along(results)) {
    entry <- results[[j]]
    if (length(entry) <= 1){
      entry_1 <- rep(NA, 32)
      entry_2 <- rep(NA, 32)
    } else {
      entry <- with(
        entry,
        c(
          halflife_nuc,
          halflife_cyt,
          deg_nuc,
          deg_cyt,
          quantiles_halflife_nuc,
          quantiles_halflife_cyt,
          quantiles_taunu,
          quantiles_lambda,
          within_quant,
          mean_coverage_nuc,
          mean_coverage_cyt,
          cv_nuc,
          cv_cyt,
          rsquared_nuc,
          rsquared_cyt
        )
      )
    }
    
    results_table <- rbind(results_table, entry)
  }
  
  column_names <- colnames(results_table)
  results_table <- data.frame(results_table)
  colnames(results_table) <- column_names
  rownames(results_table) <- genes
  print(results_table)
  
  # Create Half-Life table with reliability criteria ---------------------------
  
  elvl_regression <- expr_level_regression(
    args$input,
    time_series = args$timepoints,
    norm = "avg"
  )
  
  half_life_table <- estimation_data_table(
    rawdata = rawdata,
    estimation_results = results_table,
    elvl_regression = elvl_regression,
    timepoints = args$timepoints,
    min_cov = 30,
    max_elvl_slope_stringent = 0.002,
    max_elvl_slope_lessstringent = 0.0025,
    max_estdev = 0.33,
    max_ciq = 0.3,
    min_r2_stringent = 0.5,
    min_r2_lessstringent = 0.4
  )$all
  
  cat("Saving Results...\n")
  # Save the results
  
  save(
    results,
    results_table,
    half_life_table,
    file = args$output
  )
  cat("Computation Done.")
  
}

