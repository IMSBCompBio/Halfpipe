import gzip
from misc import r_to_c_hg38, r_to_c_m39
import multiprocessing as mp
import numpy as np
import pandas as pd
import pickle
import pysam
from scipy.stats import binom
import scipy.special
import scipy.optimize as op
import sys
from multiprocessing import set_start_method, get_context


def snp_corrected_rsummary(bamfile, reffile, snps_ai, outfile, refname_to_chr=r_to_c_hg38,
                           readfilter='singleend'):
    """
    :param bamfile:         input BAM file containing aligned reads
    :param reffile:         fasta file containing all reference sequences; an index of the same file needs to be located
                            in the same directory, having the same filename with the additional extension '.fai'
    :param snpfile:         gzipped vcf file containing SNP positions of the reference
    :param outfile:         output file to store read summary to

    :param aifile:          REDIportal file containing A-to-I editing positions of the reference (default: None)
    :param editfile:        editing sites file storing editing sites as 'ChrNumber_ChrPos' (default: None)
    :param refname_to_chr:  hash storing which reference sequence name refers to which chromosome number (needed to
                            associate the SNP sites / A-to-I sites to the reference sequences) (default: default hash)
    :param pairedend:       If True, runs function in paired-end mode. Currently, it will use Read1 reads only

    :return: void

    ....................................................................................................................

    This function collects read information for all mapped reads within a given BAM file and writes it to an output
    file.
    The information retrieved are: reference sequence name mapped to, read length, positions of insertions (position
    before insertion event, w.r.t. reference, 1-based), positions of deletions (position before deletion event, w.r.t.
    reference, 1-based), confusion table listing the number of nucleotide pairings as well as insertions and deletions
    within the read.
    Confusion table counts (nucleotide pairing counts) are corrected for SNP sites and, optionally, for editing sites.
    """

    r_summary = read_summary(bamfile, reffile, snps_ai, refname_to_chr=refname_to_chr, readfilter=readfilter)
    write_r_summary(r_summary, outfile)
    return()

def antisense_sense_rsummary(bamfile, rsummary_infile, bedfile, antisense_outfile, sense_outfile,
                             coverage=1, refname_to_chr=r_to_c_hg38, readfilter='singleend', read_orientation="reverse"):
    """
    :param bamfile:             input BAM file containing aligned reads
    :param rsummary_infile:     file storing the BAM file's read summary
    :param bedfile:             BED file containing genomic regions for which reads are to be filtered and separated
                                into antisense and sense genomic region reads
    :param antisense_outfile:   output file to store antisense genomic regions' reads summary to
    :param sense_outfile:       output file to store sense genomic regions' reads summary to

    :param coverage:            minimum coverage a genomic region has to satisfy (default: 1)
    :param refname_to_chr:      hash storing which reference sequence name refers to which chromosome number (needed to
                                associate BED chromosome names to the reference sequences) (default: default hash)

    :return: void

    ....................................................................................................................

    This function filters all reads from a BAM file that map to user-defined genomic regions and separates them into
    reads originating from antisense genomic regions and sense genomic regions. Based on the BAM file and its original
    read summary as well as a BED file specifying genomic regions, a read summary for both the antisense and sense
    genomic regions' reads is written.
    """

    genomic_regions = filter_regions(bamfile, bedfile, refname_to_chr=refname_to_chr, readfilter=readfilter, read_orientation=read_orientation)
    print("filter genomic regions", file=sys.stderr)
    genomic_regions = filter_coverage(genomic_regions, coverage=coverage)
    print("filter coverage", file=sys.stderr)
    original_rsummary = load_r_summary(rsummary_infile)
    print("load original r summary", file=sys.stderr)
    antisense_sense_separation = antisense_sense_transcripts(original_rsummary, genomic_regions)
    print("separating antisense and sense transcripts", file=sys.stderr)
    write_r_summary(antisense_sense_separation[0], antisense_outfile)
    print("writing antisense file", file=sys.stderr)
    write_r_summary(antisense_sense_separation[1], sense_outfile)
    print("writing sense file", file=sys.stderr)
    return()

def read_summary(bamfile, reffile, snps_ai, refname_to_chr=r_to_c_hg38, readfilter='singleend'):
    """
    :param bamfile: input BAM file containing aligned reads
    :param reffile: fasta file containing all reference sequences; an index of the same file needs to be located in the
                    same directory, having the same filename with the additional extension '.fai'

    :param snpfile:         gzipped vcf file containing SNP positions of the reference (default: None)
    :param aifile:          REDIportal file containing A-to-I editing positions of the reference (default: None)
    :param editfile:        editing sites file storing editing sites as 'ChrNumber_ChrPos' (default: None)
    :param refname_to_chr:  hash storing which reference sequence name refers to which chromosome number (needed to
                            associate the SNP sites / A-to-I sites to the reference sequences) (default: default hash)

    :return: a hash assigning each read's ID a list containing (in the following order):
                reference sequence name mapped to
                read length
                list of positions of insertions (position after insertion event, 1-based)
                list of positions of deletions (position of deletion event, 1-based)
                confusion table (hash) assigning a nucleotide pairing its abundance within the read's mapping

    ....................................................................................................................

    This function collects read information for all mapped reads within a given BAM file. It is designed for the
    analysis of metabolically labeled RNA transcripts.
    The information retrieved are: reference sequence name mapped to, read length, positions of insertions (position
    before insertion event, w.r.t. reference, 1-based), positions of deletions (position before deletion event, w.r.t.
    reference, 1-based), confusion table listing the number of nucleotide pairings as well as insertions and deletions
    within the read.

    The confusion table is a hash assigning each nucleotide pairing its abundance within the read's mapping. The first
    nucleotide refers to the reference, the second to the read. Pairings with 'N' are also contained (i.e., all 'NX'
    and 'XN' for some standard nucleotide X). 'XI' and 'XD' refer to insertions and deletions, respectively, where X is
    the next reference nucleotide after the insertion event or, in the case of deletions, the first nucleotide in the
    reference deleted. The counts refer to single nucleotides inserted or deleted rather than to insertion or deletion
    events.

    Optionally, confusion table counts (nucleotide pairing counts) can be corrected for SNPs and A-to-I editing sites by
    providing corresponding vcf files. Here, only single nucleotide exchanges are considered. All reference positions
    being reported to contain a SNP or being an A-to-I editing site are excluded from counting (are not comprised in
    confusion table). Also, position-wise SNP rates can be computed from an input read table file, and a maximum editing
    rate (SNP rate) can be defined for a nucleotide position to be considered.

    Read information is returned as hash, assigning a read's ID a list with its information.

    Edit 2021: For paired end, we will only consider Read 1 for the UTR conversion counts, since those reads have a
    defined priming site (at the poly-A tail). This will help us with the peak calling later on
    """

    # ..................................................................................................................
    # initializing
    # ..................................................................................................................

    bam_file = pysam.AlignmentFile(bamfile, "rb")                               # opening BAM for retrieving read IDs
    read_information = dict.fromkeys([r.query_name for r in bam_file.fetch()])  # initialize hash keys with read IDs
    bam_file.close()                                                            # closing BAM (end of file)

    bam_file = pysam.AlignmentFile(bamfile, "rb")                               # opening BAM for retrieving refseqs
    refnames = [ref["SN"] for ref in bam_file.header["SQ"]]                     # getting all reference sequence names
    bam_file.close()                                                            # closing BAM file again

    reference = pysam.FastaFile(reffile)                                        # opening reference fasta file
    if snps_ai is not None:
        editsites = pickle.load(open(snps_ai, 'rb'))                                  # load editing sites from pickle
    else:
        editsites = None
    # ..................................................................................................................
    # retrieving read information per reference sequence, since SNP and A-to-I editing sites' hash may need too much
    # memory for all reference sequences at once
    # ..................................................................................................................

    for (i, current_ref) in enumerate(refnames):

        # initializing -------------------------------------------------------------------------------------------------
        if current_ref not in refname_to_chr:
            continue

        bam_file = pysam.AlignmentFile(bamfile, "rb")           # re-opening BAM file
        if refname_to_chr and (current_ref in refname_to_chr):  # reference sequence chromosome number (is set to the
            current_chr = refname_to_chr[current_ref]           # reference sequence name if no number is given by
            if editsites is not None:
                snps_ai = {key:0 for key in editsites.keys() if key.startswith(f"{current_chr}_")}
            else:
                snps_ai = {}
        else:                                                   # 'refname_to_chr')
            current_chr = current_ref
            snps_ai = dict()

        for read in bam_file.fetch(current_ref):

            # getting read and reference metadata

            if readfilter == 'pseudosingleend': # only consider Read1 for conversion count
                if not read.is_read1:
                    continue

            ref_name = read.reference_name              # name of reference sequence mapped to read

            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue                                # (skipping unmapped and ambiguously mapped reads)
            if ref_name != current_ref:
                continue        # (skipping reads not mapped to current reference sequence)

            read_id = read.query_name                   # read ID
            read_seq = read.query_alignment_sequence    # read sequence (without hard/soft-clipped parts)
            read_length = read.query_alignment_length   # read length (without hard/soft-clipped parts)
            read_cigar = read.cigartuples               # read cigar string (as tuples)

            ref_start = read.reference_start            # reference alignment start, 0-based, including
            ref_end = read.reference_end                # reference alignment end, 0-based, excluding
            ref_seq = reference.fetch(reference=ref_name, start=ref_start, end=ref_end)     # alignment ref. seq.
            ref_seq = ref_seq.upper()                   # converting soft-masked lower-case nucleotides to upper case

            # computing and storing confusion table

            read_cigar_info = process_cigar_ntpairs(ref_seq, read_seq, read_cigar, current_chr, ref_start,
                                                    snps_ai)
            insertions = read_cigar_info[0]
            deletions = read_cigar_info[1]
            conf_table = read_cigar_info[2]

            read_information[read_id] = [ref_name, read_length, insertions, deletions, conf_table]

        bam_file.close()                                # closing BAM file (end of file)

    # ..................................................................................................................
    # returning hash storing all read information as read_id -> list of read information
    # ..................................................................................................................

    return(read_information)

def write_r_summary(r_summary, outfile):
    """
    :param r_summary:   read summary; a hash assigning each read's ID a list containing (in the following order):
                            reference sequence name mapped to
                            read orientation
                            read length
                            list of positions of T>C conversions observed (1-based)
                            list of positions of insertions (position after insertion event, 1-based)
                            list of positions of deletions (position of deletion event, 1-based)
                            hash assigning a nucleotide pairing its abundance within the read's mapping
    :param outfile:     output file to store read summary to

    :return:            void

    ....................................................................................................................

    This function writes the content of a read summary hash to a file. The file is formatted as described in the
    following:

    #read_id\n
    name_of_reference_sequence_read_is_mapped_to\n
    read_orientation\n
    read_length\n
    tab_delimited_positions_of_TC_conversions\n
    tab_delimited_positions_of_insertions\n
    tab_delimited_positions_of_deletions\n
    tab_delimited_key_value_pairs_for_nucleotide_pairing_abundances\n

    In case that there are no T>C conversions, insertions or deletions, respectively, the corresponding lines will be
    empty. The key-value pairs are written as key:value (nucleotide_pairing:abundance).
    """

    out_file = open(outfile, "w")  # opening output file
    for read_id in r_summary:  # iterating through read IDs in read summary hash
        read_record = r_summary[read_id]        # getting read ID's summary record
        out_file.write("#" + read_id + "\n")    # writing read ID
        out_file.write(read_record[0] + "\n")   # writing reference sequence mapped t
        out_file.write(str(read_record[1]) + "\n")  # writing read length
        out_file.write("\t".join([str(p) for p in read_record[2]]) + "\n")  # writing insertion positions
        out_file.write("\t".join([str(p) for p in read_record[3]]) + "\n")  # writing deletion positions
        out_file.write("\t".join([ntp + ":" + str(read_record[4][ntp]) for ntp in read_record[4]]) + "\n")  # writing
        #   nucleotide pairing abundances
    out_file.close()    # closing output file
    return ()           # returning

def filter_regions(bamfile, bedfile, refname_to_chr=r_to_c_hg38, readfilter='singleend', read_orientation="reverse"):
    """
    :param bamfile:         input BAM file containing aligned reads
    :param bedfile:         BED file containing genomic regions for which reads are to be filtered
    :param refname_to_chr:  hash storing which reference sequence name refers to which chromosome number (needed to
                            associate BED chromosome names to the reference sequences) (default: default hash)

    :return: a hash containing all genomic regions and the read IDs overlapping with these regions:
                hash keys:      genomic regions' names
                hash values:    lists containing the region's chromosome number, starting position (0-based, including),
                                ending position (0-based, excluding), strand orientation, and read IDs overlapping

    This function filters reads from a BAM file for regions as defined in a BED file, respecting read mapping
    orientation. Forwards mapped reads are only assigned to antisense (-) regions and vice versa (this is a consequence
    of the 3' sequencing approach)
    """

    # ..................................................................................................................
    # opening files and initializing
    # ..................................................................................................................

    bam_file = pysam.AlignmentFile(bamfile, "rb")       # opening BAM file with aligned reads
    bed_file = open(bedfile)                            # opening BED file with genomic regions defined
    reads_per_region = {}               # initializing hash with genomic regions and reads filtered
                                        #   (region_name -> [chr_number, region_start, region_end, strand, read_ids])

    # reference names are chromosomes: need to find genomic regions in chromosomes

    refnames = [ref["SN"] for ref in bam_file.header["SQ"]]         # getting all reference sequence names
    for i, r in enumerate(refnames):                                # replacing reference sequence names with
        if r in refname_to_chr:
            refnames[i] = refname_to_chr[r]     # chromosome numbers, if given

    genomic_regions = {}                # initializing hash with genomic regions' nucleotide positions
    for r in refnames:                  #   (chr_number -> {strand -> {nt_position -> region_name}})
        genomic_regions[r] = {"+": {}, "-": {}}

    # ..................................................................................................................
    # storing genomic regions (single nucleotide positions) in nested hash
    # chr_number -> {strand -> {nt_position -> region_name}}
    # ..................................................................................................................

    for line in bed_file:                   # iterating though BED file entries

        fields = line.strip().split("\t")                   # getting BED file entry's single fields
        chr_number = fields[0][3:]                          # genomic region chromosome number
        if chr_number not in genomic_regions:
            continue      # if chromosome number is unknown, continue

        startpos = int(fields[1])                           # genomic region starting position, 0-based, including
        endpos = int(fields[2])                             # genomic region ending position, 0-based, excluding
        name = fields[3]                                    # genomic region name
        strand = fields[5]                                  # genomic region strand orientation
        all_pos = range(startpos, endpos)                   # all positions of genomic region

        for p in all_pos:                                                   # storing the genomic region's nucleotide
            genomic_regions[chr_number][strand][p] = name                   #   positions to the genomic regions' hash

        reads_per_region[name] = [chr_number, startpos, endpos, strand]     # initialize genomic region information to
                                                                            #  reads-per-region hash

    # ..................................................................................................................
    # filtering reads from BAM file
    # region_name -> [chr_number, region_start, region_end, strand, read_ids]
    # ..................................................................................................................

    for read in bam_file.fetch():   # iterating through reads

        if readfilter == 'pseudosingleend': # we only consider Read1, since those reads map to UTRs. Read2 maps somewhere in the fragment.
            if not read.is_read1:
                continue

        read_id = read.query_name           # read ID
        ref_name = read.reference_name      # name of reference sequence mapped to read
        if ref_name in refname_to_chr:      # replacing reference name with chromosome number, if given
            ref_name = refname_to_chr[ref_name]

        if read_orientation == "reverse":
            read_o = "+" if read.is_reverse else "-" # read mapping orientation (reversed on purpose, Quant-seq 3' REV)
        elif read_orientation == "forward":
            read_o = "-" if read.is_reverse else "+" # no need to reverse oif FWD kits are used
        else:
            raise ValueError("Read orientation must be specified with 'reverse' or 'forward'")

        ref_start = read.reference_start            # reference alignment start, 0-based, including
        ref_end = read.reference_end                # reference alignment end, 0-based, excluding
        ref_positions = range(ref_start, ref_end)   # reference alignment nucleotide positions

        for p in ref_positions:                     # checking if read alignment overlaps with any genomic region
            if p in genomic_regions[ref_name][read_o]:              # overlap was found:
                gregion_name = genomic_regions[ref_name][read_o][p]     # getting genomic region name
                reads_per_region[gregion_name].append(read_id)          # storing read ID to that genomic region's list
                break                                                   # stop iteration; continue with next read

    # ..................................................................................................................
    # returning hash with genomic regions and filtered reads
    # ..................................................................................................................

    return(reads_per_region)

def filter_coverage(reads_per_region, coverage=100):
    """
    :param reads_per_region:    a hash storing genomic regions and read IDs overlapping with these regions:
                                    hash keys:      genomic regions' names
                                    hash values:    lists containing the region's chromosome number, starting position
                                                    (0-based, including), ending position (0-based, excluding), strand
                                                    orientation, and read IDs overlapping
    :param coverage:            minimum coverage a genomic region has to satisfy (default: 100)

    :return:    a hash containing only these genomic regions that hit the minimum coverage

    This function selects genomic regions according to a minimum coverage a region has to satisfy. Coverage here is
    defined as the number of reads assigned to a region.
    """

    return_hash = {}                # initializing return hash
    cov_adjusted = coverage + 4     # minimum length a genomic region's hash entry has to have

    for gr in reads_per_region:         # iterating through genomic regions

        if len(reads_per_region[gr]) >= cov_adjusted:    # checking if region hits minimum coverage
            return_hash[gr] = reads_per_region[gr]       # storing region's input hash entry to output return list

    return(return_hash)             # returning hash with selected genomic regions

def load_r_summary(infile, gene_selection=None):
    """
    :param infile:              file storing read summary
    :param gene_selection:      file storing new-line separated list of genes to load only (default: None)

    :return:        read summary; a hash assigning each read's ID a list containing (in the following order):
                        reference sequence name mapped to
                        read orientation (if stored in the input file)
                        read length
                        list of positions of insertions (position after insertion event, 1-based)
                        list of positions of deletions (position of deletion event, 1-based)
                        hash assigning a nucleotide pairing its abundance within the read's mapping

    ....................................................................................................................

    This function loads a read summary from a file. The file format should be as described in the following:

    #read_id\n
    name_of_reference_sequence_read_is_mapped_to\n
    read_orientation\n                                  (this line is optional)
    read_length\n
    tab_delimited_positions_of_TC_conversions\n
    tab_delimited_positions_of_insertions\n
    tab_delimited_positions_of_deletions\n
    tab_delimited_key_value_pairs_for_nucleotide_pairing_abundances\n

    In case that there are no T>C conversions, insertions or deletions, respectively, the corresponding lines should be
    empty. The key-value pairs should be written as key:value (nucleotide_pairing:abundance).
    """

    in_file = open(infile)  # opening read summary file
    r_summary = {}          # initialize empty read summary hash
    genes = {}              # initialize empty hash with gene names to load
    if gene_selection:      # getting gene names to load
        in_genes = open(gene_selection)
        for line in in_genes:
            genes[line.strip()] = None

    for line in in_file:    # iterating though read summary file
        read_id = line.strip()[1:]                      # getting read ID
        ref_name = next(in_file).strip()                # getting reference sequence name mapped to
        if gene_selection and ref_name not in genes:    # skipping all genes not to be loaded
            next(in_file)
            next(in_file)
            next(in_file)
            next(in_file)
            continue
        read_len = int(next(in_file).strip())           # getting read length
        ins_pos = [int(p) for p in next(in_file).strip().split("\t") if p]  # getting list of insertion positions
        del_pos = [int(p) for p in next(in_file).strip().split("\t") if p]  # getting list of deletion positions
        ntp_abundance = {}                                      # initialize hash of nucleotide pair abundances
        for key_value in next(in_file).strip().split("\t"):     # iterate through key-value pairs
            key_value = key_value.split(":")                        # splitting key-value pair
            ntp_abundance[key_value[0]] = int(key_value[1])         # storing nucleotide pair abundance

        # storing read ID's record to read_summary hash
        r_summary[read_id] = [ref_name, read_len, ins_pos, del_pos, ntp_abundance]

    in_file.close()     # closing read summary file
    return(r_summary)   # returning read summary hash

def antisense_sense_transcripts(r_summary, reads_per_region):
    """
    :param r_summary:           read summary; a hash assigning each read's ID a list containing (in the following
                                order):
                                    reference sequence name mapped to
                                    read length
                                    list of positions of T>C conversions observed (1-based)
                                    list of positions of insertions (position after insertion event, 1-based)
                                    list of positions of deletions (position of deletion event, 1-based)
                                    hash assigning a nucleotide pairing its abundance within the read's mapping
    :param reads_per_region:    a hash storing genomic regions and the read IDs overlapping with these regions:
                                    hash keys:      genomic regions' names
                                    hash values:    lists containing the region's chromosome number, starting position
                                                    (0-based, including), ending position (0-based, excluding), strand
                                                    orientation, and read IDs overlapping

    :return:    a list containing three read summaries storing antisense, sense and unspecified orientation genomic
                regions' reads
                *** NOTE *** reference sequences in the read summaries are replaced by the name of the corresponding
                genomic region

    This function separates reads from antisense (-) and sense (+) genomic regions. Reads are filtered from a read
    summary, using a hash which assigns read IDs to single genomic regions. For both antisense and sense genomic
    regions, a new read summary is returned.
    """

    # initializing read summary hashes

    ref_antisense = {}
    ref_sense = {}
    ref_unspecified = {}

    # iterating through genomic regions, separating reads by orientation

    for region in reads_per_region:

        region_entry = reads_per_region[region]     # getting genomic region's entry
        region_strand = region_entry[3]             # getting genomic region's strand orientation
        region_reads = region_entry[4:]             # getting genomic region's reads

        if (region_strand == "."):              # CASE: genomic region strand orientation not specified
            for read_id in region_reads:                # iterating through reads
                read_entry = r_summary[read_id]             # getting read's summary entry
                read_entry[0] = region                      # replacing reference sequence with genomic region
                ref_unspecified[read_id] = read_entry       # read to unspecified read summary

        elif (region_strand == "-"):            # CASE: genomic region strand orientation is antisense
            for read_id in region_reads:                # iterating through reads
                read_entry = r_summary[read_id]                 # getting read's summary entry
                read_entry[0] = region                          # replacing reference sequence with genomic region
                ref_antisense[read_id] = r_summary[read_id]     # adding all reads to forward read summary

        elif (region_strand == "+"):            # CASE: genomic region strand orientation is sense
            for read_id in region_reads:                # iterating through reads
                read_entry = r_summary[read_id]             # getting read's summary entry
                read_entry[0] = region                      # replacing reference sequence with genomic region
                ref_sense[read_id] = r_summary[read_id]     # adding all reads to reverse read summary

    # returning

    return([ref_antisense, ref_sense, ref_unspecified])

def load_snps(snpfile, chr_selected):
    """
    :param snpfile:         gzipped vcf file containing SNP positions of the reference
    :param chr_selected:    chromosome number for which SNPs are to be loaded

    :return: a hash storing single-nucleotide exchange SNPs as 'ChrNumber_ChrPos'

    This function loads single-nucleotide exchange SNPs from a gzipped vcf file for a selected chromosome number.
    """

    snps = dict()   # hash storing ChrNumber_ChrPos of SNP sites

    snp_file = gzip.open(snpfile, "rt")     # opening SNP file
    while next(snp_file).startswith("##"):  # skipping meta information (and header with last iteration)
        pass

    for line in snp_file:                   # iterating through SNP table

        fields = line.strip().split("\t")       # getting single vcf fields
        ref_chr = fields[0]                     # reference chromosome number
        if ref_chr != chr_selected:             # (skipping SNPs not from the current reference sequence)
            continue

        ref_pos = int(fields[1])                # reference chromosome position (1-based)
        ref_nt = fields[3]                      # reference nucleotides
        snp_nt = fields[4]                      # read nucleotides

        if (len(snp_nt) > 1) or (len(ref_nt) > 1):      # only considering single nucleotide exchanges
            continue
        else:
            snps[ref_chr + "_" + str(ref_pos)] = 0  # storing SNP site

    return(snps)    # returning snp hash

def load_ai_sites(aifile, chr_selected):
    """
    :param aifile:          REDIportal file containing A-to-I editing positions of the reference
    :param chr_selected:    chromosome number for which A-to-I sites are to be loaded

    :return: a hash storing A-to-I sites as 'ChrNumber_ChrPos'

    This function loads A-to-I sites from a REDIportal file for a selected chromosome number.
    """

    ai_file = open(aifile, "r")     # opening file storing A-to-I sites
    ai_hash = {}                    # hash storing ChrNumber_ChrPos of SNP sites

    for line in ai_file:            # iterating through A-to-I file
        fields = line.strip().split("\t")           # getting single A-to-I site fields
        ref_chr = fields[0][3:]                     # reference sequence name
        if ref_chr != chr_selected:
            continue        # (skipping AIs not from the selected reference sequence)

        ref_pos = fields[1]                         # reference chromosome position (1-based)
        ai_hash[ref_chr + "_" + str(ref_pos)] = 0   # storing A-to-I site

    return(ai_hash)     # returning A-to-I sites hash

def load_edit_sites(editsites_file, chr_selected):
    """
    :param editsites_file:  file storing editing sites as 'ChrNumber_ChrPos\n'
    :param chr_selected:    chromosome number for which editing sites are to be loaded

    :return: a hash storing editing sites as 'ChrNumber_ChrPos'

    This function loads editing sites from a editing sites file for a selected chromosome number.
    """

    edit_sites_file = open(editsites_file, "r")     # opening editing sites file
    edit_hash = {}                                  # hash storing editing sites as 'ChrNumber_ChrPos'

    for line in edit_sites_file:                    # iterating through editing sites file
        if line.split("_")[0] == chr_selected:          # checking if editing site belongs to selected chromosome
            edit_hash[line.strip()] = 0                     # if so, storing editing site to return hash

    edit_sites_file.close()     # closing input editing sites file
    return(edit_hash)           # returning

def process_cigar_ntpairs(ref_seq, read_seq, cigar, chr, ref_start, snp_ai_hash):
    """
    :param ref_seq:         references' nucleotide sequence the read is aligned to
    :param read_seq:        read nucleotide sequence
    :param cigar:           read cigar string
    :param chr:             reference chromosome the read is aligned to
    :param ref_start:       reference position of read alignment start
    :param snp_ai_hash:     hash storing potential SNP and A-to-I positions as ChrNumber_ChrPos

    :return:    a list containing (in the following order):
                    list of positions of insertions (position after insertion event, 1-based)
                    list of positions of deletions (position of deletion event, 1-based)
                    confusion table (hash) assigning each nucleotide pairing its abundance within the read's mapping

    This function extracts read SNP information from a given cigar string. It records positions of insertions (position
    before insertion event, w.r.t. reference, 1-based), positions of deletions (position before deletion event, w.r.t.
    reference, 1-based) and nucleotide pairing counts.
    Nucleotide pairings counts are stored in a hash, where the first nucleotide refers to the reference, the second to
    the read. 'XI' and 'XD' refer to insertions and deletions, respectively, where X is the next reference nucleotide
    after the insertion event or, in the case of deletions, the first nucleotide in the reference deleted. The counts
    refer to single nucleotides inserted or deleted rather than to insertion or deletion events.
    """

    insertions = list()     # initialize list storing insertion positions (positions w.r.t. reference,
                            # position before insertion event, 1-based)
    deletions = list()      # initialize list storing deletion positions (positions w.r.t. reference,
                            # position before deletion event, 1-based)
    conf_table = {"AA": 0, "AT": 0, "AG": 0, "AC": 0, "AN": 0, "AI": 0, "AD": 0,
                  "TA": 0, "TT": 0, "TG": 0, "TC": 0, "TN": 0, "TI": 0, "TD": 0,
                  "GA": 0, "GT": 0, "GG": 0, "GC": 0, "GN": 0, "GI": 0, "GD": 0,
                  "CA": 0, "CT": 0, "CG": 0, "CC": 0, "CN": 0, "CI": 0, "CD": 0,    # initialize confusion table
                  "NA": 0, "NT": 0, "NG": 0, "NC": 0, "NN": 0, "NI": 0, "ND": 0,
                  "RA": 0, "RT": 0, "RG": 0, "RC": 0, "RN": 0, "RI": 0, "RD": 0,# added in 2023
                  "YA": 0, "YT": 0, "YG": 0, "YC": 0, "YN": 0, "YI": 0, "YD": 0,# added in 2021    # (all counts are 0)
                  "SA": 0, "ST": 0, "SG": 0, "SC": 0, "SN": 0, "SI": 0, "SD": 0,
                  "WA": 0, "WT": 0, "WG": 0, "WC": 0, "WN": 0, "WI": 0, "WD": 0,
                  "KA": 0, "KT": 0, "KG": 0, "KC": 0, "KN": 0, "KI": 0, "KD": 0,
                  "MA": 0, "MT": 0, "MG": 0, "MC": 0, "MN": 0, "MI": 0, "MD": 0,
                  "BA": 0, "BT": 0, "BG": 0, "BC": 0, "BN": 0, "BI": 0, "BD": 0,
                  "DA": 0, "DT": 0, "DG": 0, "DC": 0, "DN": 0, "DI": 0, "DD": 0,
                  "HA": 0, "HT": 0, "HG": 0, "HC": 0, "HN": 0, "HI": 0, "HD": 0,
                  "VA": 0, "VT": 0, "VG": 0, "VC": 0, "VN": 0, "VI": 0, "VD": 0}

    # computing confusion table

    read_pos = 0        # current read position (0-based)
    ref_pos = 0         # current reference position (0-based)
    for c in cigar:     # filling confusion table with the cigar string

        operation = c[0]        # type of operation is encoded by numbers (m, ins, del)
        op_length = c[1]        # length (stretch) of operation
        snp_key_i = chr + "_"   # SNP and A-to-I sites key (initialized, incomplete)

        if operation == 4:
            continue     # skip soft-clipped positions

        if operation == 1:              # CASE: insertion
            ref_nt = ref_seq[ref_pos]                   # reference nucleotide
            table_key = ref_nt + "I"                    # generate confusion table key
            conf_table[table_key] += op_length          # incrementing insertion counter
            insertions.extend([ref_pos] * op_length)    # storing insertion positions
            read_pos += op_length                       # update read position (ref position does not change)

        elif operation == 2:            # CASE: deletion
            ref_nt = ref_seq[ref_pos]                   # reference nucleotide
            table_key = ref_nt + "D"                    # generate confusion table key
            conf_table[table_key] += op_length          # incrementing deletion counter
            deletions.extend([ref_pos] * op_length)     # storing deletion positions
            ref_pos += op_length                        # update reference position (read pos. does not change)

        else:                           # CASE: match (real match or mismatch)
            for readp, refp in zip(range(read_pos, read_pos + op_length), range(ref_pos, ref_pos + op_length)):
                if (snp_key_i + str(ref_start + 1 + refp)) in snp_ai_hash:
                    continue                            # skipping SNP / A-to-I sites
                if ref_seq[refp] == read_seq[readp]:    # real match
                    table_key = ref_seq[refp] * 2           # generate confusion table key
                    conf_table[table_key] += 1              # incrementing match counter
                else:                                   # mismatch
                    ref_nt = ref_seq[refp]                  # reference nucleotide
                    table_key = ref_nt + read_seq[readp]    # generate confusion table key
                    conf_table[table_key] += 1              # incrementing mismatch counter

            read_pos += op_length   # update read position
            ref_pos += op_length    # update reference position

    return([insertions, deletions, conf_table])     # returning

def new_convcount_table(mode, input_files, outfile):
    """
    :param mode:            type of data provided for computing conversion counts; choose 'table' for providing a
                            read table, choose 'summary' for providing a read summary
    :param input_files:     file storing the read table (mode == 'table') or the antisense and sense read summary files
                            (mode == 'summary'), respectively
    :param outfile:         output file to store conversion counts table to

    :return: void

    ....................................................................................................................

    This function computes a conversion count table from a given read table file and writes it to an output file.
    """
    cc_table = conversion_counts_table(mode, input_files)
    write_convcount_table(cc_table, outfile)
    return()

def conversion_counts_table(mode, input_files):
    """
    :param mode:            type of data provided for computing conversion counts; choose 'table' for providing a
                            read table, choose 'summary' for providing a read summary
    :param input_files:     list containing the file storing the read table (mode == 'table') or the antisense and
                            sense read summary files (mode == 'summary'), respectively

    :return:    a hash assigning each genomic region's name another hash, assigning each read ID within the genomic
                region a list containing the number of potential conversion positions and the number of actual
                conversions

    This function computes the number of potential conversion positions as well as the number of actual conversions of
    the genomic regions' reads stored within the read table file. The counts are stored as nested hash, assigning each
    genomic region's name another hash, assigning each read ID within the genomic region a list containing its counts.
    """

    region_convcounts = {}  # hash storing region_name -> {read_id -> [potential conversion positions, conversions]

    # computing conversion counts from read table

    if mode == "table":

        for gr_entry in rtable_file_generator(input_files[0]):  # iterating through read table file's genomic region
                                                                #  entries
            gr = list(gr_entry.keys())[0]                           # getting genomic region's name
            gr_counts = conversion_counts(gr_entry[gr][1])          # computing genomic region's read conversion counts
            region_convcounts[gr] = gr_counts                       # storing genomic region's conversion counts

    # computing conversion counts from read summary

    elif mode == "summary":

        antisense_summary = load_r_summary(input_files[0])      # loading antisense read summary
        sense_summary = load_r_summary(input_files[1])          # loading sense read summary

        for read in antisense_summary:        # initializing hash storing conversion counts with genomic region names
            region_convcounts[antisense_summary[read][0]] = {}
        for read in sense_summary:
            region_convcounts[sense_summary[read][0]] = {}

        for r in antisense_summary:           # iterating through antisense summary, counting conversions

            read_region = antisense_summary[r][0]               # getting read's region
            read_ntpairs = antisense_summary[r][4]              # getting read's nucleotide pairing table
            conv = read_ntpairs["AG"]                           # conversions
            pot = conv + read_ntpairs["AA"] + read_ntpairs["AT"] + read_ntpairs["AC"]   # potential conv. pos.
            region_convcounts[read_region][r] = [pot, conv]     # storing conversions

        for r in sense_summary:             # iterating through reverse summary, counting conversions

            read_region = sense_summary[r][0]                   # getting read's region
            read_ntpairs = sense_summary[r][4]                  # getting read's nucleotide pairing table
            conv = read_ntpairs["TC"]                           # conversions
            pot = conv + read_ntpairs["TA"] + read_ntpairs["TT"] + read_ntpairs["TG"]   # potential conv. pos.
            region_convcounts[read_region][r] = [pot, conv]     # storing conversions

    return(region_convcounts)       # returning

def write_convcount_table(convcount_table, outfile):
    """
    :param convcount_table:     a hash assigning each genomic region's name another hash, assigning each read ID within
                                the genomic region a list containing the number of potential conversion positions and
                                the number of actual conversions
    :param outfile:             output file to store conversion counts table to

    :return: void

    ....................................................................................................................

    This function writes the content of a conversion counts table to a file. The file is formatted as described in the
    following:

    #genomic_region\n
    >read_id\n
    number_of_potential_conversion_positions\tnumber_of_conversions\n
    """

    out_file = open(outfile, "w")       # opening output file
    for gr in convcount_table:          # iterating through genomic regions in conversion counts table
        out_file.write("#" + gr + "\n")     # writing genomic region
        gr_entry = convcount_table[gr]      # getting genomic region's entry
        for read_id in gr_entry:            # iterating through genomic region's reads
            out_file.write(">" + read_id + "\n")    # writing read ID
            out_file.write(str(gr_entry[read_id][0]) + "\t" + str(gr_entry[read_id][1]) + "\n")     # writing counts

    out_file.close()    # closing output file
    return()            # returning

def rtable_file_generator(rtable_file):
    """
    :param rtable_file:     read table file
    :return:                generator; returning one genomic region's entry in each iteration
    """

    in_file = open(rtable_file, "r")        # opening read table
    line = in_file.readline()               # reading in first line
    gr = line.strip().split("#")[1]         # initialize current genomic region (first genomic region in file)
    line = in_file.readline().strip()       # reading in second line
    gr_entry = {gr: []}                     # initialize current genomic region's entry (empty list)
    table = []                              # initialize current table (empty list)

    while line:                             # iterating through read table file until end
        if line.startswith("#"):                # next genomic region encountered
            gr_entry[gr].append(table)          # storing current conversion table to current genomic region's entry
            yield(gr_entry)                     # returning current genomic region entry
            gr = line.split("#")[1]             # setting next genomic region
            gr_entry = {gr: []}                 # initializing next genomic region's entry
            table = []                          # re-setting table
        elif line.startswith(">SNP"):       # SNP table encountered
            pass                                # nothing to do, everything was prepared in the step before
        elif line.startswith(">CONV"):      # conversion table encountered
            gr_entry[gr].append(table)          # storing current SNP table to current genomic region's entry
            table = []                          # re-setting table
        else:                               # table row
            entries = line.split("\t")          # retrieving single row entries
            # re-transforming entries into integers, Nones and strings
            entries = [entries[0]] + [int(i) if i == '0' or i == '1' else None if i == 'None' else i for i in entries[1:]]
            table.append(entries)               # storing row entries to current table
        line = in_file.readline().strip()   # loading in next line

    gr_entry[gr].append(table)              # finally, adding last conversion table to last genomic region
    yield(gr_entry)                         # finally, returning last genomic region's entry

def conversion_counts(conversion_table):
    """
    :param conversion_table:    a nested list storing read conversions:
                                    row 1: potential reference sequence conversion position, other rows: reads,
                                    column 1: read names, other columns: potential reference sequence conversion
                                        positions,
                                    entries: None if the position is not covered, masked by a deletion or SNP
                                        correction, 0 if read has no conversion at that position (may have any other
                                        SNP, though), 1 otherwise

    :return:    a hash assigning each read ID a list containing the total number of potential conversion positions and
                the number of conversions observed

    This function counts the total number of potential conversion positions as well as the number of actual conversions
    of the reads stored within a conversion table.
    """

    read_counts = {}    # hash storing read_id -> [number of potential conversion positions, number of conversions]

    for read_record in conversion_table[1:]:    # skipping first row (only contains reference genome positions)
        read_id = read_record[0]                # getting read ID
        # number of potential conversion positions (all entries unequal None and excluding read name)
        total_counts = sum([1 for i in read_record if i is not None]) - 1
        # number of conversions (all entries being 1)
        conv_counts = sum([1 for i in read_record if i == 1])
        # storing counts to hash
        read_counts[read_id] = [total_counts, conv_counts]

    return(read_counts)     # returning reads' counts

def read_summary_statistics(r_summary):
    """
    :param r_summary:   a hash assigning each read's ID a list containing information on that read (as returned by the
                        function 'read_summary'):
                            reference sequence name mapped to
                            read length
                            list of positions of insertions (position before insertion event, 1-based)
                            list of positions of deletions (position before deletion event, 1-based)
                            hash assigning a nucleotide pairing its abundance within the read's mapping

    :return: a list containing (in the following order):
                total number of reads
                hash: fragment length -> ratio of reads of that length
                relative position of insertions
                relative position of deletions
                a hash storing for each nucleotide pairing:
                    number of pairings -> % reads with that many pairings
                    % reference nucleotides paired with this read nucleotide (only considering SNPs, not indels)
                    % reads that have at least one such pairing and additionally contain at least one insertion
                    % reads that have at least one such pairing and additionally contain at least one deletion
                    relative position of additional insertions
                    relative position of additional deletions
                    hash: read length -> % reads of that length (of reads with at least one such pairing)

    ....................................................................................................................

    This function computes a summary statistics on read information collected for all mapped reads within a given BAM
    file as returned by the function 'read_summary'. Just as 'read_summary', it is designed for the analysis of
    metabolically labeled RNA transcripts.

    The statistics includes:

    total number of reads
    fragment length distribution as hash: fragment length -> % reads of that length
    relative position of insertions
    relative position of deletions
    a hash storing for each nucleotide pairing:
        number of pairings distribution as hash: number of pairings -> % reads with that many pairings
        % reference nucleotides being paired with the corresponding read nucleotide (only considering SNPs, not indels)
        % reads that have at least one such pairing and additionally contain at least one insertion
        % reads that have at least one such pairing and additionally contain at least one deletion
        relative position of additional insertions
        relative position of additional deletions
        fragment length distr. (of reads with at least one such pairing) as hash: read length -> % reads of that length

    The percentage of reference nucleotides being paired with the corresponding read nucleotide are computed from the
    total number of reference nucleotides being paired to any other nucleotide, i.e. insertion and deletion events are
    not taken into account here.
    Relative positions of insertions and deletions are computed as the average relative position within an read and
    then averaged over all reads. Position of an inserted or deleted nucleotide is the position of the last matched
    nucleotide before the insertion or deletion event. Each nucleotide is considered individually rather than taking
    into account only the insertion or deletion event; therefore, larger indels have higher impact on an read's average
    relative indel position.
    """

    # ..................................................................................................................
    # initializing
    # ..................................................................................................................

    total_reads = len(r_summary)            # total number of reads
    read_lengths = {}                       # hash storing read length -> % reads of that length
    ins_positions = 0                       # relative position of insertions
    del_positions = 0                       # relative position of deletions

    # stats to be recorded for each nucleotide pairing:
        # number of pairings -> % reads with that many pairings,
        # % reference nucleotides being paired with the corresponding read nucleotide (doesn't consider indels),
        # % reads that additionally contain insertions, % reads that additionally contain deletions,
        # average relative position of additional insertion, average relative position of additional deletion,
        # read length -> % reads of that length
    ntpair_stats = {
        "AA": [{0: 0}, 0, 0, 0, 0, 0, {}], "AT": [{0: 0}, 0, 0, 0, 0, 0, {}], "AG": [{0: 0}, 0, 0, 0, 0, 0, {}],
        "AC": [{0: 0}, 0, 0, 0, 0, 0, {}], "AI": [{0: 0}, 0, 0, 0, 0, 0, {}], "AD": [{0: 0}, 0, 0, 0, 0, 0, {}],
        "TA": [{0: 0}, 0, 0, 0, 0, 0, {}], "TT": [{0: 0}, 0, 0, 0, 0, 0, {}], "TG": [{0: 0}, 0, 0, 0, 0, 0, {}],
        "TC": [{0: 0}, 0, 0, 0, 0, 0, {}], "TI": [{0: 0}, 0, 0, 0, 0, 0, {}], "TD": [{0: 0}, 0, 0, 0, 0, 0, {}],
        "GA": [{0: 0}, 0, 0, 0, 0, 0, {}], "GT": [{0: 0}, 0, 0, 0, 0, 0, {}], "GG": [{0: 0}, 0, 0, 0, 0, 0, {}],
        "GC": [{0: 0}, 0, 0, 0, 0, 0, {}], "GI": [{0: 0}, 0, 0, 0, 0, 0, {}], "GD": [{0: 0}, 0, 0, 0, 0, 0, {}],
        "CA": [{0: 0}, 0, 0, 0, 0, 0, {}], "CT": [{0: 0}, 0, 0, 0, 0, 0, {}], "CG": [{0: 0}, 0, 0, 0, 0, 0, {}],
        "CC": [{0: 0}, 0, 0, 0, 0, 0, {}], "CI": [{0: 0}, 0, 0, 0, 0, 0, {}], "CD": [{0: 0}, 0, 0, 0, 0, 0, {}]}
    nt_pairs = ["AA", "AT", "AG", "AC", "AI", "AD", "TA", "TT", "TG", "TC", "TI", "TD", "GA", "GT", "GG", "GC",
                "GI", "GD", "CA", "CT", "CG", "CC", "CI", "CD"]     # keys for nucleotide pairing statistics hash

    # ..................................................................................................................
    # iterating through read summary, filling read summary statistics
    # ..................................................................................................................

    ins = 0     # total number of insertion affected reads (needed for averaging)
    de = 0      # total number of deletion affected reads (needed for averaging)

    for read_id in r_summary:    # iterating through reads

        read_records = r_summary[read_id]    # all information on read
        read_len = read_records[1]                  # length of read
        read_ins = read_records[2]                  # read insertion positions
        read_del = read_records[3]                  # read deletion positions
        read_conf = read_records[4]                 # read confusion table

        # filling read length distribution -----------------------------------------------------------------------------


        if read_len in read_lengths:
            read_lengths[read_len] += 1
        else:
            read_lengths[read_len] = 1

        # counting up insertions and relative insertion positions ------------------------------------------------------

        is_ins = 0          # variable indicating whether insertions are present
        pos_ins = 0         # initialize average relative position of read's insertions
        if read_ins:
            is_ins = 1                                              # tagging that insertions are present in read
            pos_ins = sum(read_ins) / (len(read_ins) * read_len)    # compute average relative pos. of read's insertions

        ins += is_ins               # incrementing insertion reads counter
        ins_positions += pos_ins    # adding average relative position to insertion relative position counter

        # counting up deletions and relative deletion positions --------------------------------------------------------

        is_del = 0          # variable indicating whether insertions are present
        pos_del = 0         # initialize average relative position of read's deletions
        if read_del:
            is_del = 1                                              # tagging that deletions are present in read
            pos_del = sum(read_del) / (len(read_del) * read_len)    # compute average relative pos. of read's deletions

        de += is_del                # incrementing deletion reads counter
        del_positions += pos_del    # adding average relative position to deletion relative position counter

        # filling nucleotide pairing statistics ------------------------------------------------------------------------

        for pair in nt_pairs:       # iterating through nucleotide pairings
            pair_frequency = read_conf[pair]                # getting nt pairing's frequency from confusion table

            # counting up nucleotide pairing frequencies
            if pair_frequency in ntpair_stats[pair][0]:     # if existent: incrementing pairing frequency counter
                ntpair_stats[pair][0][pair_frequency] += 1
            else:                                           # if not yet existent: initialize frequency counter with 1
                ntpair_stats[pair][0][pair_frequency] = 1

            # (following statistics refer to observing at least one nt pairing)
            if pair_frequency > 0:

            # counting up additional insertions and relative insertion positions
                ntpair_stats[pair][2] += is_ins
                ntpair_stats[pair][4] += pos_ins

            # counting up additional deletions and relative deletion position
                ntpair_stats[pair][3] += is_ins
                ntpair_stats[pair][5] += pos_del

            # filling read length distribution
                if read_len in ntpair_stats[pair][6]:
                    ntpair_stats[pair][6][read_len] += 1
                else:
                    ntpair_stats[pair][6][read_len] = 1

    # ..................................................................................................................
    # computing averages and ratios (yet, absolute numbers stored)
    # ..................................................................................................................

    # counting number of reference A, T, C, G nucleotides for % nucleotide pairings
    ref_A = sum([sum([f * ntpair_stats[p][0][f] for f in ntpair_stats[p][0]]) for p in ["AA", "AT", "AC", "AG"]])
    ref_T = sum([sum([f * ntpair_stats[p][0][f] for f in ntpair_stats[p][0]]) for p in ["TA", "TT", "TC", "TG"]])
    ref_C = sum([sum([f * ntpair_stats[p][0][f] for f in ntpair_stats[p][0]]) for p in ["CA", "CT", "CC", "CG"]])
    ref_G = sum([sum([f * ntpair_stats[p][0][f] for f in ntpair_stats[p][0]]) for p in ["GA", "GT", "GC", "GG"]])

    # % reads for fragment length distribution
    for i in read_lengths:
        read_lengths[i] = read_lengths[i] / total_reads
    # average relative insertion position
    if ins:
        ins_positions /= ins
    # average relative deletion position
    if de:
        del_positions /= de

    # nucleotide pairing information
    for pair in ntpair_stats:

        # % reference nucleotides paired with corresponding read nucleotide (not considering indels)
        if pair[1] != "I" and pair[1] != "D":
            ref_N = ref_A                                                       # reference nucleotide (A)
            if pair[0] == "T":
                ref_N = ref_T                                    # reference nucleotide (T)
            elif pair[0] == "C":
                ref_N = ref_C                                  # reference nucleotide (C)
            elif pair[0] == "G":
                ref_N = ref_G                                  # reference nucleotide (G)
            pair_counts = sum(
                [f * ntpair_stats[pair][0][f] for f in ntpair_stats[pair][0]])  # nt pair counts
            ntpair_stats[pair][1] = pair_counts / ref_N                         # % ref nt paired with read nt

        # % reads for nt pair frequency distribution
        pair_observed = total_reads - ntpair_stats[pair][0][0]  # number of reads with at least one such nt pair
        if not pair_observed:                                   # if there is no such nt pair observed,
            ntpair_stats[pair][0][0] /= total_reads                 # only adapt 0 - occurrence information
            continue                                                # and continue
        for f in ntpair_stats[pair][0]:
            ntpair_stats[pair][0][f] /= total_reads

        # % reads with additional insertion
        ntpair_stats[pair][2] /= pair_observed
        # % reads with additional deletions
        ntpair_stats[pair][3] /= pair_observed
        # average relative insertion position
        if ntpair_stats[pair][2]:
            ntpair_stats[pair][4] /= ntpair_stats[pair][2]
        # average relative deletion position
        if ntpair_stats[pair][3]:
            ntpair_stats[pair][5] /= ntpair_stats[pair][3]
        # % reads for read length distribution
        for i in ntpair_stats[pair][6]:
            ntpair_stats[pair][6][i] /= pair_observed

    # ..................................................................................................................
    # returning
    # ..................................................................................................................

    return([total_reads, read_lengths, ins_positions, del_positions, ntpair_stats])

def new_conveff_table(cctable_file, outfile, common_regions, base_error, sequencing_error, coverage=200, iterations=100,
                      n_labelingsites=30):
    """
    :param cctable_file:        input conversion counts table file
    :param outfile:             output file to store conversion efficiency table to
    :param base_error:          probability of observing a conversion due to a sequencing error
    :param sequencing_error:    probability that the converted nucleotide is masked by a sequencing error

    :param fix_efficiency:      optional; a value to fix for conversion efficiency (as a consequence, only new/total
                                ratio will be estimated) (default: None
    :param new_ini:             optional; initial value for new/total ratio (default: None)
    :param coverage:            minimum coverage a genomic region has to satisfy for estimation (default: 50)
    :param iterations:          number of iterations to use for efficiency estimation (default: 100)

    :return: void

    ....................................................................................................................

    This function computes a conversion efficiency table from a given conversion counts table file and writes it to an
    output file. In detail, the ratio of newly synthesized transcripts and the conversion efficiency are estimated for
    each genomic region using an EM algorithm. Optionally, a minimum coverage threshold can be set a genomic region has
    to fulfill to be considered for estimation.
    """
    conveff_table = conversion_efficiency_table(cctable_file, common_regions, base_error, sequencing_error,
                                                n_labelingsites=n_labelingsites, coverage=coverage, iterations=iterations)
    write_conveff_table(conveff_table, outfile)
    return()

def conversion_efficiency_table(cctable_file, common_regions, base_error, sequencing_error, n_labelingsites=30, coverage=100, iterations=100):
    """
    :param cctable_file:        file storing genomic regions' read conversion counts table
    :param base_error:          probability of observing a conversion due to a sequencing error
    :param sequencing_error:    probability that the converted nucleotide is masked by a sequencing error

    :param fix_efficiency:  optional; a value to fix for conversion efficiency (as a consequence, only new/total ratio
                            will be estimated) (default: None
    :param new_ini:         optional; initial value for new/total ratio (default: None)
    :param coverage:        optional; minimum coverage a genomic region has to satisfy for estimation (default: 50)
    :param iterations:      optional; number of iterations to use for efficiency estimation (default: 100)

    :return:    a hash assigning each genomic region's name a list containing two sublists (the final estimations being
                the last sublist entries, respectively):
                    the first sublist contains the sequence of estimated newly synthesized by total RNA ratios
                    the second sublist contains the sequence of estimated conversion efficiencies

    This function estimates the ratio of newly synthesized transcripts as well as the conversion efficiency for each
    genomic region, based on a given conversion counts table file. Optionally, a minimum coverage threshold can be set
    a genomic region has to fulfill to be considered for estimation (all regions not fulfilling the minimum coverage
    will be assigned a 'None' in the returned hash).
    """

    # loading in genomic regions' read conversion counts from file
    # (hash storing genomic_region -> {read_id -> [potential conversion positions, conversions]})
    cc_table = load_convcount_table(cctable_file)
    # hash storing genomic_region -> [[newly synthesized ratio estimations], [conversion efficiency estimations]]
    estimations = {}

    est = global_efficiency_estimation(
        cc_table,
        base_error,
        sequencing_error,
        common_regions,
        iterations=iterations,
        coverage=coverage,
        n_labelingsites=n_labelingsites
    )  # estimation

    for genomic_region in cc_table:
        estimations[genomic_region] = [est[0][genomic_region], est[1]]

    return estimations


def write_conveff_table(conveff_table, outfile):
    """
    :param conveff_table:   a hash assigning each genomic region's name a list containing two sublists:
                                the first sublist contains the sequence of estimated newly synthesized by total RNA
                                    ratios
                                the second sublist contains the sequence of estimated conversion efficiencies
    :param outfile:         output file to store conversion efficiency table to

    :return: void

    ....................................................................................................................

    This function writes the content of a conversion efficiencies table to a file. The file is formatted as described
    in the following:

    #genomic_region\n
    >new_by_total\n
    tab-separated listing of new-by-total ratio estimations\n
    >conv_eff\n
    tab-separated listing of conversion efficiency estimations\n
    """

    out_file = open(outfile, "w")  # opening output file
    for gr in conveff_table:  # iterating through genomic regions in conversion efficiency table

        out_file.write("#" + gr + "\n")  # writing genomic region
        out_file.write(">new_by_total\n")  # writing new-by-total ratios headline
        if conveff_table[gr][0] is None:
            out_file.write("None\n")
        else:
            for i in conveff_table[gr][0][:-1]:  # iterating through new-by-total entries
                out_file.write(str(i) + "\t")  # writing entry
            out_file.write(str(conveff_table[gr][0][-1]) + "\n")

        out_file.write(">conv_eff\n")  # writing conversion efficiencies headline
        if conveff_table[gr][0] is None:
            out_file.write("None\n")
        else:
            for i in conveff_table[gr][1][:-1]:  # iterating through conversion efficiency entries
                out_file.write(str(i) + "\t")  # writing entry
            out_file.write(str(conveff_table[gr][1][-1]) + "\n")

    out_file.close()  # closing output file
    return()            # returning

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

def global_efficiency_estimation(
        cc_table,
        base_error,
        sequencing_error,
        common_regions,
        iterations=50,
        coverage=200,
        n_labelingsites=30,
        fix_efficiency=None
):
    """
    :param r_conversions:       a hash assigning each read ID a list containing the total number of potential conversion
                                positions and the number of conversions observed
    :param base_error:          probability of observing a conversion due to a sequencing error
    :param sequencing_error:    probability that the converted nucleotide is masked by a sequencing error

    :param fix_efficiency:      optional; a value to fix for conversion efficiency (as a consequence, only new/total
                                ratio will be estimated) (default: None
    :param new_ini:             optional; initial value for new/total ratio (default: None)
    :param iterations:          optional; number of iterations to use for efficiency estimation (default: 100)

    :return:    a list containing two sublists: the first sublist contains the sequence of estimated newly synthesized
                by total RNA ratios, the second sublist contains the sequence of estimated conversion efficiencies (the
                final estimations being the last list entries); in case no estimation could be made (labeled read counts
                too low), the list contains two 'None' values instead

    This function estimates the conversion efficiency from a collection of reads, using an EM algorithm that in parallel
    estimates the new/total ratio as well.
    """

    # ..................................................................................................................
    # collecting all combinations: (number of potential conversion positions and number of conversions) read counts
    # hash: number potential conversion pos. -> {number conversions -> number of reads}
    # ..................................................................................................................
    newtotals = {}
    all_pos_conv_numbers = {}

    conv_1_global = 0
    conv_2_global = 0
    conv_all_global = 0
    npos_global = 0

    valid_regions = []
    for genomic_region in cc_table:

        r_conversions = cc_table[genomic_region]
        valid_utrs = np.where(np.array(list(r_conversions.values()))[:, 0] >= n_labelingsites)[0]
        if np.shape(valid_utrs)[0] < coverage or genomic_region not in common_regions:  # minimum coverage threshold is not fulfilled:
            pos_conv_numbers = {}  # initializing hash
            for r in r_conversions:  # iterating through reads
                p = r_conversions[r][0]  # number of potential conversion positions
                c = r_conversions[r][1]  # number of conversions
                if p in pos_conv_numbers:  # number of positions already stored in hash
                    if c in pos_conv_numbers[p]:  # number of conversions already stored in hash
                        pos_conv_numbers[p][c] += 1  # incrementing number of reads in existing entry
                    else:  # number of conversions not yet stored in hash
                        pos_conv_numbers[p][c] = 1  # setting new entry
                else:  # number of positions not yet stored in hash
                    pos_conv_numbers[p] = {c: 1}  # setting new entry
            all_pos_conv_numbers[genomic_region] = pos_conv_numbers
            # ..................................................................................................................
            # initial estimations of 'new' (ratio of newly synthesized RNA transcripts) and 'e' (conversion efficiency)
            # initialize 'new' with the ratio of labeled by total RNA transcripts
            # initialize 'e' with an efficiency estimation excluding base and sequencing errors, or base error
            # ..................................................................................................................

            conv_1 = 0  # number of reads with one conversion
            conv_2 = 0  # number of reads with two conversions
            conv_all = 0  # number of reads with one or more conversions (number of labeled reads)
            for p in pos_conv_numbers:  # iterating over hash storing combination counts, counting up labeled reads
                for c in pos_conv_numbers[p]:  # iterating over conversion numbers
                    if c == 1:  # conversion number 1
                        conv_1 += pos_conv_numbers[p][c]  # incrementing counter for one-time converted reads
                    elif c == 2:  # conversion number 2
                        conv_2 += pos_conv_numbers[p][c]  # incrementing counter for two-time converted reads
                    elif c != 0:
                        conv_all += pos_conv_numbers[p][c]  # incrementing counter for multiple-times converted reads
            conv_all += conv_1 + conv_2  # adding 1- and 2-converted read counts to
            # one-or-more-converted read counts

            new = conv_all / len(r_conversions)  # new: labeled/total transcript counts ratio
            newtotals[genomic_region] = [new]

        else:
            valid_regions.append(genomic_region)
            r_conversions = {key: r_conversions[key] for key in list(np.array(list(r_conversions.keys()))[valid_utrs])}

            pos_conv_numbers = {}  # initializing hash
            for r in r_conversions:  # iterating through reads
                p = r_conversions[r][0]  # number of potential conversion positions
                c = r_conversions[r][1]  # number of conversions
                if p in pos_conv_numbers:  # number of positions already stored in hash
                    if c in pos_conv_numbers[p]:  # number of conversions already stored in hash
                        pos_conv_numbers[p][c] += 1  # incrementing number of reads in existing entry
                    else:  # number of conversions not yet stored in hash
                        pos_conv_numbers[p][c] = 1  # setting new entry
                else:  # number of positions not yet stored in hash
                    pos_conv_numbers[p] = {c: 1}  # setting new entry
            all_pos_conv_numbers[genomic_region] = pos_conv_numbers
        # ..................................................................................................................
        # initial estimations of 'new' (ratio of newly synthesized RNA transcripts) and 'e' (conversion efficiency)
        # initialize 'new' with the ratio of labeled by total RNA transcripts
        # initialize 'e' with an efficiency estimation excluding base and sequencing errors, or base error
        # ..................................................................................................................

            conv_1 = 0      # number of reads with one conversion
            conv_2 = 0      # number of reads with two conversions
            conv_all = 0    # number of reads with one or more conversions (number of labeled reads)
            for p in pos_conv_numbers:      # iterating over hash storing combination counts, counting up labeled reads
                for c in pos_conv_numbers[p]:   # iterating over conversion numbers
                    if c == 1:                      # conversion number 1
                        conv_1 += pos_conv_numbers[p][c]    # incrementing counter for one-time converted reads
                    elif c == 2:                    # conversion number 2
                        conv_2 += pos_conv_numbers[p][c]    # incrementing counter for two-time converted reads
                    elif c != 0:
                        conv_all += pos_conv_numbers[p][c]  # incrementing counter for multiple-times converted reads
            conv_all += conv_1 + conv_2     # adding 1- and 2-converted read counts to
            # one-or-more-converted read counts

            new = conv_all / len(r_conversions)  # new: labeled/total transcript counts ratio
            newtotals[genomic_region] = [new]

            conv_1_global += conv_1
            conv_2_global += conv_2
            conv_all_global += conv_all

            if (conv_1 != 0 and conv_2 != 0):  # e: else, ratio of 2 and 1-times converted, excluding base errors
                npos_global += list(pos_conv_numbers.keys())[0] + 2


    # initial values for 'e' and 'new'


    e = base_error                              # e: initialize with base error
    if fix_efficiency:                          # e: if given, fixed efficiency
        e = fix_efficiency
    elif (conv_1_global != 0 and conv_2_global != 0):         # e: else, ratio of 2 and 1-times converted, excluding base errors
        ratio = conv_2_global/conv_1_global                           # ratio of 2- and 1-times converted read counts
        e = ratio / ((npos_global - 1) / 2 + ratio)            # initial 'e' estimation

    #print(e)

    # ..................................................................................................................
    # iteratively estimating 'new' and 'e' using EM algorithm
    # ..................................................................................................................

    e_estimation_sequence = [e]         # sequence of 'e' estimations during EM algorithm iterations

    genomic_regions = all_pos_conv_numbers.keys()

    print(f"{len(valid_regions)}/{len(list(genomic_regions))}  3'UTRs ({len(valid_regions)/len(list(genomic_regions)) * 100}%) used for Labeling Efficiency Estimation.", file=sys.stdout)
    for i in range(iterations):  # iteratively re-estimating

        # conversion probability at one position for newly synthesized transcripts
        a = e * (1 - sequencing_error - base_error) + base_error
        #print(a)

        # P(H|O ; theta) probability of hidden variable (pre-existing or newly synthesized transcript) given observation
        # (read number of conversions and number of potential conversion positions)
        # hash: number of potential conversion positions -> {number of conversions -> P(H | O; theta)}

        w0_global = 0  # sum over all P(H=0|O=c; theta)
        w1_global = 0  # sum over all P(H=1|O=c; theta)
        v_global = 0  # sum over all ( P(H=1|O=c; theta) * c )
        z_global = 0  # sum over all ( P(H=1|O=c; theta) * (p-c) )

        for genomic_region in genomic_regions:

            #r_conversions = cc_table[genomic_region]
            pos_conv_numbers = all_pos_conv_numbers[genomic_region]
            new = newtotals[genomic_region][-1]
            prob_h0 = dict.fromkeys(pos_conv_numbers, {})       # initialize hash for H = 0 (pre-existing)
            prob_h1 = dict.fromkeys(pos_conv_numbers, {})       # initialize hash for H = 1 (newly synthesized)

            for p in pos_conv_numbers:      # iterating over number of potential conversion positions
                for c in pos_conv_numbers[p]:   # iterating over number of conversions
                    combi = scipy.special.binom(p, c)       # combinations to draw 'c' out of 'p' (binomial coefficient)
                    h0 = combi * (base_error ** c) * ((1 - base_error) ** (p - c)) * (1 - new)   # P(O|H)*P(H) for H = 0
                    h1 = combi * (a ** c) * ((1 - a) ** (p - c)) * new                           # P(O|H)*P(H) for H = 1
                    denominator = h0 + h1                                                        # P(O); denominator
                    prob_h0[p][c] = h0 / denominator        # computing and storing P(H|O; theta) for H = 0
                    prob_h1[p][c] = h1 / denominator        # computing and storing P(H|O; theta) for H = 1

            # sums over P(H|O; theta) needed for subsequent computations

            w0 = 0      # sum over all P(H=0|O=c; theta)
            w1 = 0      # sum over all P(H=1|O=c; theta)
            v = 0       # sum over all ( P(H=1|O=c; theta) * c )
            z = 0       # sum over all ( P(H=1|O=c; theta) * (p-c) )

            for p in pos_conv_numbers:  # iterating over number of potential conversion positions
                for c in pos_conv_numbers[p]:   # iterating over number of conversions
                    counts = pos_conv_numbers[p][c]     # read counts
                    h0 = prob_h0[p][c]                  # corresponding P(H=0|O=c)
                    h1 = prob_h1[p][c]                  # corresponding P(H=1|O=c)
                    w0 += counts * h0           # incrementing w0
                    w1 += counts * h1           # incrementing w1
                    v += counts * h1 * c        # incrementing v
                    z += counts * h1 * (p - c)  # incrementing z

            if genomic_region in valid_regions:
                w0_global += w0
                w1_global += w1
                v_global += v
                z_global += z

            # re-estimating 'new' and 'e', appending to estimation sequences

            new = w1 / (w0 + w1)
            newtotals[genomic_region].append(new)

        if fix_efficiency:          # fixed conversion efficiency:
            e = fix_efficiency          # setting e to fixed value
        else:
            if (v_global + z_global) == 0:
                e = base_error
            else:
                e = 1 / (1 - base_error - sequencing_error) * (v_global / (v_global + z_global) - base_error)
        if e < 0:
            e = 0
        if e > 1:
            e = 1
        e_estimation_sequence.append(e)

    def likelihood(pi, n_conversions, n_uridines, base_error, le):
        prob = np.sum(np.log(((1 - pi) * binom.pmf(n_conversions, n_uridines, base_error)) + (
                    pi * binom.pmf(n_conversions, n_uridines, le))))
        return -prob

    newtotals = {}
    a = e * (1 - sequencing_error - base_error) + base_error
    for genomic_region in genomic_regions:
        reads = cc_table[genomic_region]
        n_uridines = np.array(list(reads.values()))[:, 0]  # vector containing potential conversion sites
        n_conversions = np.array(list(reads.values()))[:, 1]  # vector containing observed conversion sites
        mode = op.minimize(likelihood, x0=a, args=(n_conversions, n_uridines, base_error, a), bounds=((0, 1),),
                           method="L-BFGS-B").x[0]
        newtotals[genomic_region] = [mode]

    # ..................................................................................................................
    # returning
    # ..................................................................................................................

    return([newtotals, e_estimation_sequence])

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

# Wrapper --------------------------------------------------------------------------------------------------------------

def read_summary_wrapper(infile):

    print(f"summary wrapper says infile 17 is {infile[17]}", file=sys.stderr)

    if infile[17] == 'human':
        refname_to_chr = r_to_c_hg38
    elif infile[17] == 'mouse':
        refname_to_chr = r_to_c_m39
        print('mouse organism used', file=sys.stderr)
    else:
        raise ValueError('organism must be human or mouse')

    snp_corrected_rsummary(
        bamfile=infile[0],
        reffile=infile[1],
        snps_ai=infile[3],
        outfile=infile[4],
        readfilter=infile[12],
        refname_to_chr=refname_to_chr
    )
    print('snp corrected', file=sys.stderr)
    return None


def antisense_sense_wrapper(infile):
    if infile[17] == 'human':
        refname_to_chr = r_to_c_hg38
    elif infile[17] == 'mouse':
        refname_to_chr = r_to_c_m39
        print('mouse organism used', file=sys.stderr)
    else:
        raise ValueError('organism must be human or mouse')

    antisense_sense_rsummary(
        bamfile=infile[0],
        rsummary_infile=infile[4],
        bedfile=infile[2],
        antisense_outfile=infile[5],
        sense_outfile=infile[6],
        coverage=1,
        readfilter=infile[12],
        read_orientation=infile[16],
        refname_to_chr=refname_to_chr)

    print('antisense summary', file=sys.stderr)
    return None


def new_convcount_table_wrapper(infile):
    outfile = infile[7]
    infiles = infile[5:7]
    new_convcount_table('summary', infiles, outfile)
    return None


def base_and_sequencing_error(infile):
    if infile[8] == 'control' or str(infile[8]) == '0':
        stats = read_summary_statistics(load_r_summary(infile[5]))
        base_error_antisense = stats[4]['AG'][1]
        sequencing_error_antisense = stats[4]['GA'][1] + stats[4]['GT'][1] + stats[4]['GC'][1]
        stats = read_summary_statistics(load_r_summary(infile[6]))
        base_error_sense = stats[4]['TC'][1]
        sequencing_error_sense = stats[4]['CA'][1] + stats[4]['CT'][1] + stats[4]['CG'][1]
        base_error = np.mean([base_error_sense, base_error_antisense])
        print(f"Base Error: {round(base_error*100, 2)}%", file=sys.stderr)
        sequencing_error = np.mean([sequencing_error_sense, sequencing_error_antisense])
        print(f"Sequencing Error: {round(sequencing_error * 100, 2)}%", file=sys.stderr)
        pickle.dump([base_error, sequencing_error], open(infile[9], 'wb'))
        return None
    else:
        return None

def get_commonlyexpressed_regions(conv_files):
    common_regions = []
    cctables = []
    for cctable_file in conv_files:
        cc_table = load_convcount_table(cctable_file)
        cctables.append(cc_table)
        common_regions.append(list(cc_table.keys()))
    common_regions = list(set.intersection(*map(set, common_regions)))

    return common_regions

def new_conveff_table_wrapper(infile):
    controlfile = infile[9]
    baseerror, sequencingerror = pickle.load(open(controlfile, 'rb'))

    new_conveff_table(
        infile[7],
        infile[11],
        infile[16],
        baseerror,
        sequencingerror,
        coverage=infile[18],
        iterations=100,
        n_labelingsites=infile[17]
    )
    return None

def editsitecorrection_stats(infile):
    """
    This function generates a CSV-file that contains information about all conventional base conversion statistics.

    :param infile: Internal pyhton list for parallelisation.
    :return: None, CSV-file gets directly saved in the ratioestimation subfolder.
    """
    antisensestats_fullcorrection = read_summary_statistics(load_r_summary(infile[5]))
    sensestats_fullcorrection = read_summary_statistics(load_r_summary(infile[6]))

    convs = [
        'AA', 'AC', 'AG', 'AT',
        'GG', 'GA', 'GC', 'GT',
        'CC', 'CA', 'CG', 'CT',
        'TT', 'TA', 'TC', 'TG'
    ]

    stats = [
        [antisensestats_fullcorrection[4][conv][1] for conv in convs],
        [sensestats_fullcorrection[4][conv][1] for conv in convs]
    ]

    stats = pd.DataFrame(stats)
    stats.columns = convs
    stats['name'] = [
        'antisense_fullcorrection',
        'sense_fullcorrection'
    ]

    stats.to_csv(infile[20], index=False)
    ax = stats.drop(['AA', 'TT', 'CC', 'GG'], axis=1).plot(x='name',
                    kind='bar',
                    stacked=False,
                    title='Nucleotide Conversion Rates')
    ax.set_xticklabels(stats['name'], rotation=45)
    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(infile[21], dpi=300)
    return None

def EstimateRatios(samples, control, config):

    if config['params']['filterstrategy'] == 'pairedend' or config['params']['filterstrategy'] == 'pseudosingleend':
        indicator1 = 2
        indicator2 = 3
    else:
        indicator1 = 1
        indicator2 = 2


    input_files = [
        [
            f"{config['output']}filtering/{sample[indicator2]}_mapped_filtered.bam",
            config['input']['refgenome'],
            config['input']['bed'],
            f"{config['output']}readpreprocess/{control}_alleditsites.pkl",
            f"{config['output']}ratioestimation/{sample[indicator2]}_rsummary.txt",
            f"{config['output']}ratioestimation/{sample[indicator2]}_antisensersummary.txt",
            f"{config['output']}ratioestimation/{sample[indicator2]}_sensersummary.txt",
            f"{config['output']}ratioestimation/{sample[indicator2]}_convcount.txt",
            sample[indicator1],
            f"{config['output']}ratioestimation/{control}_readsummarystats.pkl",
            f"{config['output']}ratioestimation/{sample[indicator2]}_conveff.txt",
            f"{config['output']}ratioestimation/{sample[indicator2]}_convratio.txt",
            config['params']['filterstrategy'],
            f"{config['output']}ratioestimation/{sample[indicator2]}_conversionstats.csv",
            f"{config['output']}ratioestimation/{sample[indicator2]}_conversionstats_barplot.png",
            f"{config['output']}ratioestimation/{sample[indicator2]}_labelingefficiencies.pkl",
            config['params']['readorientation'],
            config['params']['organism']

        ]
        for sample in samples
    ]



    print("Creating Conversion Count Tables...", file=sys.stderr)
    with mp.get_context("spawn").Pool(processes=int(config['params']['cores'])) as pool:
        pool.map(read_summary_wrapper, input_files)

    print("Separating Reads By Strand Orientation...", file=sys.stderr)
    with mp.get_context("spawn").Pool(processes=int(config['params']['cores'])) as pool:
        pool.map(antisense_sense_wrapper, input_files) # problem is here, fsr it stops computing ?????

    print("Creating Conversion Count Tables...", file=sys.stderr)
    with mp.get_context("spawn").Pool(processes=int(config['params']['cores'])) as pool:
        pool.map(new_convcount_table_wrapper, input_files)

    conv_files = []
    for sample in samples:
        conv_files.append(f"{config['output']}ratioestimation/{sample[indicator2]}_convcount.txt")
    common_regions = get_commonlyexpressed_regions(conv_files)

    input_files = [
        [
            f"{config['output']}filtering/{sample[indicator2]}_mapped_filtered.bam",
            config['input']['refgenome'],
            config['input']['bed'],
            f"{config['output']}readpreprocess/{control}_alleditsites.pkl",
            f"{config['output']}ratioestimation/{sample[indicator2]}_rsummary.txt",
            f"{config['output']}ratioestimation/{sample[indicator2]}_antisensersummary.txt",
            f"{config['output']}ratioestimation/{sample[indicator2]}_sensersummary.txt",
            f"{config['output']}ratioestimation/{sample[indicator2]}_convcount.txt",
            sample[indicator1],
            f"{config['output']}ratioestimation/{control}_readsummarystats.pkl",
            f"{config['output']}ratioestimation/{sample[indicator2]}_conveff.txt",
            f"{config['output']}ratioestimation/{sample[indicator2]}_convratio.txt",
            config['params']['filterstrategy'],
            f"{config['output']}ratioestimation/{sample[indicator2]}_conversionstats.csv",
            f"{config['output']}ratioestimation/{sample[indicator2]}_conversionstats_barplot.png",
            f"{config['output']}ratioestimation/{sample[indicator2]}_labelingefficiencies.pkl",
            common_regions,
            config['emparams']['uridinethreshold'],
            config['emparams']['minimumcoverage'],
            config['params']['organism']

        ]
        for sample in samples
    ]


    print("Estimating Labeling Efficiencies And New/Total Ratios...", file=sys.stderr)
    with mp.Pool(processes=int(config['params']['cores'])) as pool:
        pool.map(base_and_sequencing_error, input_files)

    print("Estimation of Nuclear and Cytosolic New/Total Ratios...", file=sys.stdout)
    with mp.Pool(processes=int(config['params']['cores'])) as pool:
        pool.map(new_conveff_table_wrapper, input_files)

    print("\n", file=sys.stdout)

    return None

