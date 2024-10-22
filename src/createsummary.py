from misc import r_to_c_hg38
import numpy as np
import pysam

def new_summary_table(outfile, description, libsize, cctable_file, cetable_file=None, coverage=1,
                      refname_to_chr=r_to_c_hg38):
    """
    :param outfile:             output file to write summary table to
    :param description:         a string describing the data processed in the summary matrix
    :param libsize:             library size of the data processed in the summary matrix
    :param cctable_file:        conversion counts table file

    :param cetable_file:        conversion efficiency table file (default: None)
    :param coverage:            minimum coverage a genomic region has to satisfy (default: 1)
    :param refname_to_chr:      hash storing which reference sequence name refers to which chromosome number (needed to
                                associate BED chromosome names to the reference sequences) (default: default hash)

    :return:    a hash assigning each genomic region's name a list containing (in the following order) the description,
                unlabeled read counts, labeled read counts

    ....................................................................................................................

    This function computed a summary table for a given a conversion counts table and conversion efficiency table file
    and writes the table to an output file. If no conversion efficiency table file is given, estimated conversion
    efficiencies and newly synthesized ratios are set to -1. Genomic regions to be included in the matrix can be
    selected by defining a coverage minimum (as defined by parameter 'coverage'). Coverage here is defined as the
    number of reads associated to a genomic region.

    The output file is of tab-delimited format:
    region_name\tdescription\tlibrary_size\ttotal_reads\tlabeled_reads\taverage_potential_conversion_positions\t
    conv_efficiency\tnewly_synthesized_ratio\n
    """

    stable = summary_table(description, libsize, cctable_file, cetable_file=cetable_file,
                           coverage=coverage, refname_to_chr=refname_to_chr)    # computing summary table
    write_summary_table(stable, outfile)                                        # writing summary table
    return()

def summary_table(description, libsize, cctable_file, cetable_file=None, coverage=1, refname_to_chr=r_to_c_hg38):
    """
    :param description:         a string describing the data processed in the summary matrix
    :param libsize:             library size of the data processed in the summary matrix
    :param cctable_file:        conversion counts table file

    :param cetable_file:        conversion efficiency table file (default: None)
    :param coverage:            minimum coverage a genomic region has to satisfy (default: 1)
    :param refname_to_chr:      hash storing which reference sequence name refers to which chromosome number (needed to
                                associate BED chromosome names to the reference sequences) (default: default hash)

    :return:    a nested hash, the outer hash assigning each genomic region's name an inner hash assigning the
                description a list containing library size, total read counts, labeled read counts, conversion
                efficiency estimation and newly synthesized transcripts ratio estimation:
                region_name -> description -> [libsize, total, labeled, average number of potential conversion positions,
                                               conv. efficiency, newly synthesized ratio]

    ....................................................................................................................

    This function records, for each genomic region, the number of total and labeled reads as well as the estimated
    conversion efficiencies and newly synthesized transcripts ratios, given a conversion counts table file and a
    conversion efficiency table file. If no conversion efficiency table file is given, estimated conversion efficiencies
    and newly synthesized ratios are set to -1. A description of the data set processed as well as the total library
    size has to be given, too.

    Optionally, only such regions can be chosen for evaluation that fulfill a certain coverage minimum (as defined by
    parameter 'coverage'). Coverage here is defined as the number of reads associated to a genomic region.
    """

    stable = {}                                         # initializing summary table (empty hash)
    cc_table = load_convcount_table(cctable_file)       # loading in conversion counts table
    if cetable_file:
        ce_table = load_conveff_table(cetable_file)     # if given, loading in conversion efficiency table
    else:
        ce_table = None                                 # else, setting ce_table to None

    # computing read counts per region

    for gr in cc_table:             # iterating through genomic regions

        gr_cc_entry = cc_table[gr]          # getting region's cc table entry
        gr_total = len(cc_table[gr])        # getting genomic region's total read counts
        if gr_total < coverage:             # skip genomic regions not fulfilling minimum coverage
            continue

        gr_convpos = 0                      # total number of potential conversion positions in all reads of the region
        gr_labeled = 0                      # number of labeled reads
        for read_id in gr_cc_entry:         # iterating through reads
            gr_convpos += gr_cc_entry[read_id][0]   # counting up total number of potential conversion positions
            if gr_cc_entry[read_id][1] > 0:         # counting up reads with at least one conversion
                gr_labeled += 1
        gr_cp_avg = gr_convpos / gr_total   # average total number of potential conversion positions

        if ce_table:                        # if given, getting conversion efficiency and newly synthesized ratio
            gr_ce_entry = ce_table[gr]          # getting region's ce table entry
            if gr_ce_entry[1] is not None:          # if successfully computed, getting region's estimations
                gr_conveff = gr_ce_entry[1][-1]     # getting region's estimated conversion efficiency
                gr_ratio = gr_ce_entry[0][-1]       # getting region's estimated newly synthesized transcripts ratio
            else:                               # else, setting estimations to -1
                gr_conveff = -1
                gr_ratio = -1
        else:                               # else, setting efficiency and ratio to -1
            gr_conveff = -1                     # setting region's estimated conversion efficiency to -1
            gr_ratio = -1                       # setting region's estimated newly synthesized transcripts ratio to -1

        # storing genomic region's read counts to summary table hash
        stable[gr] = {description: [libsize, gr_total, gr_labeled, gr_cp_avg, gr_conveff, gr_ratio]}

    # returning

    return(stable)

def write_summary_table(stable, outfile):
    """
    :param stable:      a nested hash, the outer hash assigning each genomic region's name an inner hash assigning the
                        description a list containing library size, total read counts, labeled read counts, conversion
                        efficiency estimation and newly synthesized transcripts ratio estimation:
                        region_name -> description -> [libsize, total, labeled, average potential conversion positions,
                                                       conv. efficiency, newly ratio]
    :param outfile:     output file to write summary matrix to

    :return: void

    This function writes a summary table to a file (according to the following format):
    region_name\tdescription\tlibrary_size\ttotal_reads\tlabeled_reads\taverage_potential_conversion_positions\t
    conv_efficiency\tnewly_synthesized_ratio\n
    """

    out_file = open(outfile, "w")   # opening output file

    for gr in stable:               # iterating through summary table hash, writing entries to output file
        for d in stable[gr]:
            entry = stable[gr][d]
            out_file.write(gr + "\t" + d + "\t" + "\t".join([str(e) for e in entry]) + "\n")

    out_file.close()        # closing output file
    return ()               # returning

def load_convcount_table(infile):
    """
    :param infile:  file storing conversion counts table

    :return:    a hash assigning each genomic region's name another hash, assigning each read ID within the genomic
                region a list containing the number of potential conversion positions and the number of actual
                conversions

    ....................................................................................................................

    This function loads a conversion counts table from a file. The file format should be as described in the following:

    #genomic_region\n
    >read_id\n
    number_of_potential_conversion_positions\tnumber_of_conversions\n
    """

    in_file = open(infile)  # opening conversion counts table file
    convcounts_table = {}   # initialize empty conversion counts table hash

    current_gr = next(in_file).strip()[1:]      # storing current genomic region being processed
    gr_entry = {}                               # initialize current genomic region's entry
    current_read = next(in_file).strip()[1:]    # storing current read ID being processed
    for line in in_file:                        # iterating through conversion counts table file
        if line.startswith("#"):                    # new genomic region encountered
            convcounts_table[current_gr] = gr_entry     # storing former genomic region
            current_gr = line.strip()[1:]               # setting new genomic region
            gr_entry = {}                               # re-setting genomic region's entry
        elif line.startswith(">"):                  # new read ID encountered
            current_read = line.strip()[1:]             # setting new read ID
        else:                                       # counts encountered
            gr_entry[current_read] = [int(c) for c in line.strip().split("\t")]     # storing counts
    convcounts_table[current_gr] = gr_entry     # storing last genomic region, too

    in_file.close()             # closing input file
    return(convcounts_table)    # returning

def load_conveff_table(infile):
    """
    :param infile:  file storing conversion efficiency table

    :return:    a hash assigning each genomic region's name a list containing two sublists:
                    the first sublist contains the sequence of estimated newly synthesized by total RNA ratios
                    the second sublist contains the sequence of estimated conversion efficiencies

    ....................................................................................................................

    This function loads a conversion efficiency table from a file. The file format should be as described in the
    following:

    #genomic_region\n
    >new_by_total\n
    tab-separated listing of new-by-total ratio estimations\n
    >conv_eff\n
    tab-separated listing of conversion efficiency estimations\n
    """

    in_file = open(infile)  # opening conversion efficiency table file
    conveff_table = {}      # initialize empty conversion efficiency table hash

    current_gr = next(in_file).strip()[1:]      # storing current genomic region being processed
    gr_entry = []                               # initialize current genomic region's entry
    for line in in_file:                        # iterating through conversion efficiency table file
        if line.startswith("#"):                    # new genomic region encountered
            conveff_table[current_gr] = gr_entry        # storing former genomic region
            current_gr = line.strip()[1:]               # setting new genomic region
            gr_entry = []                               # re-setting genomic region's entry
        elif line.startswith(">"):                  # headline for new-by-total ratios or conv. eff. encountered
            pass                                        # nothing to do, everything was prepared in the step before
        else:                                       # single entries encountered
            if line.startswith("None"):
                gr_entry.append(None)                                   # transforming entries to None values
            else:
                entries = [float(i) for i in line.strip().split("\t")]  # transforming entries to floats
                gr_entry.append(entries)
    conveff_table[current_gr] = gr_entry        # storing last genomic region, too

    in_file.close()         # closing input file
    return(conveff_table)   # returning

def merged_summary_table(mode, stable_files, outfile):
    """
    :param mode:            mode to be executed: choose 'integrate' for integrating tables' entries, choose 'merge'
                            for merging tables into one matrix
    :param stable_files:    list of summary table files
    :param outfile:         output file to store integrated / merged summary table to

    :return: void

    ....................................................................................................................

    This function takes multiple summary tables and either integrates or merges them into one.

    Integration will integrate all entries that are of equal genomic region name and data description. In detail,
    library size as well as read counts will be summed up, average conversion positions as well as conversion effi-
    ciencies and newly synthesized ratios will be weightedly averaged. Only genomic region name and data description
    combinations present in all input summary tables will be considered.
    Merging will collect genomic regions' entries. Only genomic regions that are present in all input summary tables
    will be collected.

    The integrated / merged summary table will be sorted by genomic regions.
    """

    # initializing

    out_file = open(outfile, "w")                                           # opening output file
    all_stables = [load_summary_table(file) for file in stable_files]       # loading all summary tables
    first_table = all_stables[0]                                            # getting first table
    shared_regions = []                 # list storing genomic regions that occur in all summary tables

    # filtering genomic regions (and data descriptions) that occur in all summary tables

    if mode == "integrate":     # CASE: integrate matrices

        for gr in first_table:              # iterating through regions (of first summary table)

            shared = True                   # storing if region + [descriptions] occurs in all tables (set to True)
            for des in first_table[gr]:     # iterating through region's data descriptions
                for tab in all_stables[1:]:     # iterating through remaining matrices
                    if gr not in tab:                       # checking if region is missing in matrix
                        shared = False                      # 'shared' set to False (doesn't occur in all matrices)
                    elif des not in tab[gr]:                # checking if description is missing in matrix
                        shared = False                      # 'shared' set to False (doesn't occur in all matrices)
            if shared:                      # if region + [descriptions] occurs in all matrices,
                for des in first_table[gr]:    # store region and descriptions to selected regions
                    shared_regions.append([gr, des])

    elif mode == "merge":       # CASE: merge matrices

        for gr in first_table:              # iterating through regions (of first summary matrix)

            shared = True                           # storing if region occurs in all matrices (set to True)
            for tab in all_stables[1:]:             # iterating through remaining matrices
                if gr not in tab:                       # checking if region is missing in matrix
                    shared = False                      # 'shared' set to False (doesn't occur in all matrices)

            if shared:                              # if region occurs in all matrices:
                shared_regions.append(gr)           # store region to selected regions

    # integrating summary tables sorted by region name

    if mode == "integrate":

        for gr_des in shared_regions:       # iterating through shared genomic regions + descriptions
            gr = gr_des[0]                      # genomic region name
            des = gr_des[1]                     # data description
            libsize = 0                         # summed library size (initialized with 0)
            total = 0                           # summed total read counts (initialized with 0)
            labeled = 0                         # summed labeled read counts (initialized with 0)
            convpos = 0                         # weighted summed potential conversion positions
            conveff = 0                         # weighted summed conversion efficiencies (initialized with 0)
            ratios = 0                          # weighted summed newly synthesized ratios (initialized with 0)

            for tab in all_stables:             # iterating through summary tables
                tab_total = tab[gr][des][1]             # table entry's total number of reads
                libsize += tab[gr][des][0]              # summing up library size
                total += tab_total                      # summing up total reads
                labeled += tab[gr][des][2]              # summing up labeled reads
                convpos += tab_total * tab[gr][des][3]  # summing up weighted potential conversion positions
                conveff += tab_total * tab[gr][des][4]  # summing up weighted conversion efficiencies
                ratios += tab_total * tab[gr][des][5]   # summing up weighted newly synthesized ratios

            convpos /= total                    # averaging potential conversion positions
            conveff /= total                    # averaging conversion efficiencies
            ratios /= total                     # averaging newly synthesized ratios

            out_file.write(
                gr + "\t" + des + "\t" +
                "\t".join([str(i) for i in [libsize, total, labeled, convpos,
                                            conveff, ratios]]) + "\n")              # writing entry

    # merging summary matrices sorted by region name

    elif mode == "merge":

        for gr in shared_regions:               # iterating through shared genomic regions
            for tab in all_stables:                 # iterating through summary matrices
                for d in tab[gr]:                       # iterating through genomic region's entries
                    entry = tab[gr][d]                                                      # getting entry
                    out_file.write(
                        gr + "\t" + d + "\t" + "\t".join([str(e) for e in entry]) + "\n")   # writing entry

    # returning

    return()

def load_summary_table(infile):
    """
    :param infile:  file storing the summary table

    :return:    a nested hash, the outer hash assigning each genomic region's name an inner hash assigning the
                description a list containing library size, total read counts, labeled read counts, conversion
                efficiency estimation and newly synthesized transcripts ratio estimation:
                region_name -> description -> [libsize, total, labeled, average potential conversion positions,
                                               conv. efficiency, newly ratio]

    This function loads a summary table from a file. The file format should be as described in the following:
    region_name\tdescription\tlibrary_size\ttotal_reads\tlabeled_reads\tconv_efficiency\tnewly_synthesized_ratio\n
    """

    stable = {}                     # initializing summary table as empty hash
    in_file = open(infile, "r")     # opening summary table file

    for line in in_file:            # iterating through summary matrix file
        entry = line.strip().split("\t")    # getting a genomic region's entries
        gr_name = entry[0]                  # getting a genomic region's name

        if gr_name in stable:  # genomic region already stored:
            stable[gr_name][entry[1]] = [float(e) for e in entry[2:]]       # adding description's entry
        else:  # genomic region not yet stored:
            stable[entry[0]] = {entry[1]: [float(e) for e in entry[2:]]}    # setting up new genomic region's entry

    in_file.close()     # closing summary matrix file
    return(stable)      # returning summary matrix

def CreateSummaryFiles(input, config):

    if config['params']['model'] == "onecompartment":
        input = input.iloc[:, 0].tolist()
    if config['params']['model'] == "twocompartment":
        input = input.iloc[:, 0].tolist() + input.iloc[:, 1].tolist()
    summaryfiles = []
    print(input)
    for file in input:
        task = pysam.view('-c', f"{config['output']}/filtering/{file}_mapped_filtered.bam")
        libsize = int(task) / 2
        summaryfile = f"{config['output']}/summaryfiles/{file}_summarytable_globalem.tsv"
        new_summary_table(summaryfile,
                          file, libsize,
                          f"{config['output']}/ratioestimation/{file}_convcount.txt",
                          cetable_file=f"{config['output']}/ratioestimation/{file}_convratio.txt")
        summaryfiles.append(summaryfile)

    merged_summary_table(mode="merge", stable_files=summaryfiles, outfile=f"{config['output']}/summaryfiles/{config['input']['summaryname']}.tsv")
    summaryfiles = np.array(summaryfiles)
    print(summaryfiles)
    summaryfiles = np.delete(summaryfiles, [0, int(np.shape(summaryfiles)[0]/2)])
    merged_summary_table(mode="merge", stable_files=summaryfiles, outfile=f"{config['output']}/summaryfiles/{config['input']['summaryname']}_nocontrol.tsv")




