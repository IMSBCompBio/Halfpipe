import gzip
from misc import r_to_c_hg38, chr_numbers_hs
import multiprocessing as mp
import os
import pickle
import pysam
import sys

def read_table(mode, bamfile, reffile, bedfile, snpfile=None, coverage=1, refname_to_chr=r_to_c_hg38, readfilter='singleend'):
    """
    :param mode:        mode to be executed; choose 'snp_all' to compute SNP tables reporting all kind of conversions,
                        choose 'snp_exclusive' to compute SNP tables reporting all but labeling-specific conversions,
                        choose 'conv' to compute conversion tables only, use 'snp_all_summary', 'snp_exclusive_summary',
                        'conv_summary' to compute summarized tables that report SNP and total read counts for each
                        position, choose 'both' or 'both_summary' to compute both a SNP (snp_exclusive) and a conversion
                        table
    :param bamfile:     input BAM file containing aligned reads
    :param reffile:     fasta file containing all reference sequences; an index of the same file needs to be located in
                        the same directory, having the same filename with the additional extension '.fai'
    :param bedfile:     BED file containing genomic regions for which reads are to be filtered

    :param snpfile:         gzipped vcf file containing SNP positions of the reference (default: None)
    :param coverage:        minimum coverage a genomic region has to satisfy (default: 1)
    :param refname_to_chr:  hash storing which reference sequence name refers to which chromosome number (needed to
                            associate the SNP sites / A-to-I sites to the reference sequences) (default: default hash)

    :return:    a hash assigning each genomic region's name a list containing a SNP and a conversion table (implemented
                as nested lists, the first index corresponding to rows and the second to columns):
                    SNP table: row 1: genomic region's reference sequence positions (0-based), other rows: reads,
                        column 1: read names, other columns: reference sequence positions, entries: None if the position
                        is not covered by the read, 0 if read has no SNP at that position, SNP nucleotide otherwise
                    conversion table: row 1: potential reference sequence conversion position, other rows: reads,
                        column 1: read names, other columns: potential reference sequence conversion positions,
                        entries: None if the position is not covered, masked by a deletion or SNP correction, 0 if read
                        has no conversion at that position (may have any other SNP, though), 1 otherwise

    This function records SNPs as well as conversions for each genomic region, resolved by reads. SNPs and conversions
    are stored in separate tables, with rows corresponding to reads and columns to potential SNP or conversion posi-
    tions, respectively, in the reference genome. The first row always contains reference sequence positions, the
    first column always contains read IDs (in case of summary tables, the read IDs are set to 'total_counts' and
    'event_counts'). The first row + first column entry is set to None. SNPs are denoted as None if the position is not
    covered by the read, 0 in case no SNP is present, and the SNP nucleotide otherwise (deletions are denoted with 'D',
    insertions with 'I' in additional columns at the end of the table). Conversions are stored as None in case that the
    position is not covered by the read, masked by a deletion or SNP correction, 0 in case no conversion is present,
    and as 1 otherwise.

    Optionally, conversion positions can be corrected for SNP positions by providing a SNP vcf file. In this case, all
    positions listed to be prone to SNPs will be reported as None in the table.
    Optionally, a minimum coverage can be specified by which genomic regions are filtered; only regions hitting this
    minimum coverage will be processed and included into the read table hash.
    """

    # ..................................................................................................................
    # initializing
    # ..................................................................................................................

    print("\nINITIALIZING")

    reference = pysam.FastaFile(reffile)  # opening reference file
    bam_file = pysam.AlignmentFile(bamfile, "rb")  # opening BAM for retrieving refseqs
    refnames = [ref["SN"] for ref in bam_file.header["SQ"]]  # getting all reference sequence names
    bam_file.close()  # closing BAM file again

    # filtering reads and genomic regions ------------------------------------------------------------------------------

    genomic_regions = filter_regions(  # filtering reads by genomic region
        bamfile, bedfile, refname_to_chr=refname_to_chr, readfilter=readfilter)
    for gr in list(genomic_regions.keys()):  # filtering genomic regions by strand orientation
        if genomic_regions[gr][3] == ".":  # (filter out unspecified strand orientations)
            del genomic_regions[gr]
    if coverage:  # filtering genomic regions by coverage (keeping regions hitting minimum coverage threshold)
        genomic_regions = filter_coverage(genomic_regions, coverage=coverage)
    read_regions = {}  # reverting genomic regions' hash: read_ID -> genomic region and strand orientation
    for gr in genomic_regions:
        for read_id in genomic_regions[gr][4:]:
            read_regions[read_id] = [gr, genomic_regions[gr][3]]

    # initialize SNP and conversion table related hashes ---------------------------------------------------------------

    # SNP table related hashes

    regions_minpos = {}  # hash storing region_name -> minimum position covered (0-based, including)
    regions_maxpos = {}  # hash storing region_name -> maximum position covered (0-based, excluding)
    region_insinfo = {}  # hash storing region_name -> {position -> {read_ids with ins at that position}}
    read_snpinfo = {}  # hash storing read_ID -> [starting position (0-based, incl), SNPs per position]

    if mode.startswith("snp") or mode.startswith("both"):
        for gr in genomic_regions:
            regions_minpos[gr] = sys.maxsize
            regions_maxpos[gr] = 0
            region_insinfo[gr] = {}

    # conversion table related hashes

    regions_convpos = {}  # hash storing region_name -> hash of potential conversion positions (0-based)
    snps = {}  # hash storing ChrNumber_ChrPos of reported SNPs
    read_convinfo = {}  # hash storing read_ID -> [starting position (0-based, incl), conv per position]

    if mode.startswith("conv") or mode.startswith("both"):
        for gr in genomic_regions:
            regions_convpos[gr] = {}

    # initializing hash storing read tables

    read_tables = dict.fromkeys([k for k in genomic_regions])

    # loading reported SNPs' hash (if needed)
    if mode.startswith("conv") or mode.startswith("both"):
        print("loading ALL SNPs ...", file=sys.stderr)
        if snpfile:
            # snps = load_snps(snpfile, current_chr)
            snps = load_snps(snpfile)  # load all SNPs at once
            print(snps)
        print(" done", file=sys.stderr)

    # ..................................................................................................................
    # collecting reads' SNP and conversion information
    # (per reference sequence, since SNP sites' hash may need too much memory for all eference sequences at once)
    # ..................................................................................................................

    for (i, current_ref) in enumerate(refnames):

        # initializing -------------------------------------------------------------------------------------------------


        if current_ref not in refname_to_chr:
            continue
        print(f"PROCESSING REFERENCE {current_ref} ...", file=sys.stderr)
        bam_file = pysam.AlignmentFile(bamfile, "rb")  # re-opening BAM file
        if refname_to_chr and (current_ref in refname_to_chr):  # reference sequence chromosome number (is set to the
            current_chr = refname_to_chr[current_ref]  # reference sequence name if no number is given by
        else:  # 'refname_to_chr')
            current_chr = current_ref

        # iterating through all reads in the BAM file getting reads' SNPs and conversions per positions-----------------

        print("computing SNP and conversion tables ...")
        for read in bam_file.fetch(current_ref):

            if readfilter == 'pseudosingleend':
                if not read.is_read1:
                    continue

            read_id = read.query_name  # read ID
            if read_id not in read_regions:
                continue  # skipping reads not mapped to any genomic region
            if read.reference_name != current_ref:
                continue  # skipping reads not mapped to current refseq
            read_region = read_regions[read_id][0]  # read genomic region
            region_strand = read_regions[read_id][1]  # read genomic region's strand (ref forward or reverse)
            if region_strand == ".":
                continue  # skipping reads with unspecified strand

            read_seq = read.query_alignment_sequence  # read sequence (without hard/soft-clipped parts)
            ref_start = read.reference_start  # reference alignment start, 0-based, including
            ref_end = read.reference_end  # reference alignment end, 0-based, excluding
            ref_seq = reference.fetch(  # reference alignment seq (soft-masked lower-case patched)
                reference=read.reference_name, start=ref_start, end=ref_end).upper()

            update_region_borders(read_region, regions_minpos, regions_maxpos,
                                  ref_start, ref_end)  # updating genomic regions' min and max pos
            (converted_nt, conversion_nt) = labeling_ntpair(region_strand)  # defining conversion (labeling) nt pair

            # getting and storing read SNP and conversion information

            construction_mode = mode.split("_summary")[0]
            read_snpconv = process_cigar_snps_conv(construction_mode, read.cigartuples, converted_nt, conversion_nt,
                                                   read_id, read_seq, current_chr, read_region, ref_seq, ref_start,
                                                   region_insinfo, regions_convpos, snps=snps)  # read SNP and conv info

            read_snpinfo[read_id] = [ref_start] + read_snpconv[0]
            read_convinfo[read_id] = [ref_start] + read_snpconv[1]

        # closing BAM file (end of file)

        bam_file.close()

    # ..................................................................................................................
    # constructing genomic regions' read tables and storing them
    # ..................................................................................................................

    print("\nCONSTRUCTING SNP AND CONVERSION TABLES ...")

    for gr in genomic_regions:  # iterating through all genomic regions, building read tables
        gr_snptable = []  # initializing SNP table for current genomic region
        gr_convtable = []  # initializing conversion table for current genomic region

        if mode.startswith("snp") or mode.startswith("both"):  # --- SNP table ---
            gr_snptable = construct_snp_table(genomic_regions[gr][4:], read_snpinfo, region_insinfo[gr],
                                              regions_minpos[gr], regions_maxpos[gr])
        if mode.startswith("conv") or mode.startswith("both"):  # --- conversion table ---
            gr_convtable = construct_conv_table(genomic_regions[gr][4:], read_convinfo,
                                                list(regions_convpos[gr].keys()))
        if mode.endswith("summary"):  # --- summary tables ---
            construction_mode = mode.split("_")[0]
            table_summaries = construct_snp_conv_summary_table(construction_mode, gr_snptable, gr_convtable)
            gr_snptable = table_summaries[0]
            gr_convtable = table_summaries[1]

        read_tables[gr] = [gr_snptable, gr_convtable]  # adding genomic region's table to return hash

    print(" done")

    # ..................................................................................................................
    # returning hash storing genomic regions' read tables
    # ..................................................................................................................

    return (read_tables)


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

    bam_file = pysam.AlignmentFile(bamfile, "rb")  # opening BAM file with aligned reads
    bed_file = open(bedfile)  # opening BED file with genomic regions defined
    reads_per_region = {}  # initializing hash with genomic regions and reads filtered
    #   (region_name -> [chr_number, region_start, region_end, strand, read_ids])

    # reference names are chromosomes: need to find genomic regions in chromosomes

    refnames = [ref["SN"] for ref in bam_file.header["SQ"]]  # getting all reference sequence names
    for i, r in enumerate(refnames):  # replacing reference sequence names with
        if r in refname_to_chr:
            refnames[i] = refname_to_chr[r]  # chromosome numbers, if given

    genomic_regions = {}  # initializing hash with genomic regions' nucleotide positions
    for r in refnames:  # (chr_number -> {strand -> {nt_position -> region_name}})
        genomic_regions[r] = {"+": {}, "-": {}}

    # ..................................................................................................................
    # storing genomic regions (single nucleotide positions) in nested hash
    # chr_number -> {strand -> {nt_position -> region_name}}
    # ..................................................................................................................

    for line in bed_file:  # iterating though BED file entries

        fields = line.strip().split("\t")  # getting BED file entry's single fields
        chr_number = fields[0][3:]  # genomic region chromosome number
        if chr_number not in genomic_regions:
            continue  # if chromosome number is unknown, continue

        startpos = int(fields[1])  # genomic region starting position, 0-based, including
        endpos = int(fields[2])  # genomic region ending position, 0-based, excluding
        name = fields[3]  # genomic region name
        strand = fields[5]  # genomic region strand orientation
        all_pos = range(startpos, endpos)  # all positions of genomic region

        for p in all_pos:  # storing the genomic region's nucleotide
            genomic_regions[chr_number][strand][p] = name  # positions to the genomic regions' hash

        reads_per_region[name] = [chr_number, startpos, endpos, strand]  # initialize genomic region information to
        #  reads-per-region hash

    # ..................................................................................................................
    # filtering reads from BAM file
    # region_name -> [chr_number, region_start, region_end, strand, read_ids]
    # ..................................................................................................................

    for read in bam_file.fetch():  # iterating through reads

        if readfilter == 'pseudosingleend':  # only consider Read1, since those reads map to UTRs
            if not read.is_read1:
                continue

        read_id = read.query_name  # read ID
        ref_name = read.reference_name  # name of reference sequence mapped to read
        if ref_name in refname_to_chr:  # replacing reference name with chromosome number, if given
            ref_name = refname_to_chr[ref_name]

        if read_orientation == "reverse":
            read_o = "+" if read.is_reverse else "-"  # read mapping orientation (reversed on purpose, Quant-seq 3' REV)
        else:
            read_o = "-" if read.is_reverse else "+" # no need to reverse if FWD kits are used

        ref_start = read.reference_start  # reference alignment start, 0-based, including
        ref_end = read.reference_end  # reference alignment end, 0-based, excluding
        ref_positions = range(ref_start, ref_end)  # reference alignment nucleotide positions

        for p in ref_positions:  # checking if read alignment overlaps with any genomic region
            if p in genomic_regions[ref_name][read_o]:  # overlap was found:
                gregion_name = genomic_regions[ref_name][read_o][p]  # getting genomic region name
                reads_per_region[gregion_name].append(read_id)  # storing read ID to that genomic region's list
                break  # stop iteration; continue with next read

    # ..................................................................................................................
    # returning hash with genomic regions and filtered reads
    # ..................................................................................................................

    return (reads_per_region)


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

    return_hash = {}  # initializing return hash
    cov_adjusted = coverage + 4  # minimum length a genomic region's hash entry has to have

    for gr in reads_per_region:  # iterating through genomic regions

        if len(reads_per_region[gr]) >= cov_adjusted:  # checking if region hits minimum coverage
            return_hash[gr] = reads_per_region[gr]  # storing region's input hash entry to output return list

    return (return_hash)  # returning hash with selected genomic regions


def load_snps(snpfile, chr_selected=None):
    """
    :param snpfile:         gzipped vcf file containing SNP positions of the reference
    :param chr_selected:    chromosome number for which SNPs are to be loaded

    :return: a hash storing single-nucleotide exchange SNPs as 'ChrNumber_ChrPos'

    This function loads single-nucleotide exchange SNPs from a gzipped vcf file for a selected chromosome number.
    """

    snps = dict()  # hash storing ChrNumber_ChrPos of SNP sites

    snp_file = gzip.open(snpfile, "rt")  # opening SNP file
    while next(snp_file).startswith("##"):  # skipping meta information (and header with last iteration)
        pass

    for line in snp_file:  # iterating through SNP table

        fields = line.strip().split("\t")  # getting single vcf fields
        ref_chr = fields[0]  # reference chromosome number

        if chr_selected is not None:
            if ref_chr != chr_selected:  # (skipping SNPs not from the current reference sequence)
                continue

        ref_pos = int(fields[1])  # reference chromosome position (1-based)
        ref_nt = fields[3]  # reference nucleotides
        snp_nt = fields[4]  # read nucleotides

        if (len(snp_nt) > 1) or (len(ref_nt) > 1):  # only considering single nucleotide exchanges
            continue
        else:
            snps[ref_chr + "_" + str(ref_pos)] = 0  # storing SNP site

    return (snps)  # returning snp hash


def load_ai_sites(aifile, chr_selected):
    """
    :param aifile:          REDIportal file containing A-to-I editing positions of the reference
    :param chr_selected:    chromosome number for which A-to-I sites are to be loaded

    :return: a hash storing A-to-I sites as 'ChrNumber_ChrPos'

    This function loads A-to-I sites from a REDIportal file for a selected chromosome number.
    """

    ai_file = open(aifile, "r")  # opening file storing A-to-I sites
    ai_hash = {}  # hash storing ChrNumber_ChrPos of SNP sites

    for line in ai_file:  # iterating through A-to-I file
        fields = line.strip().split("\t")  # getting single A-to-I site fields
        ref_chr = fields[0][3:]  # reference sequence name
        if ref_chr != chr_selected:
            continue  # (skipping AIs not from the selected reference sequence)

        ref_pos = fields[1]  # reference chromosome position (1-based)
        ai_hash[ref_chr + "_" + str(ref_pos)] = 0  # storing A-to-I site

    return (ai_hash)  # returning A-to-I sites hash


def load_edit_sites(editsites_file, chr_selected):
    """
    :param editsites_file:  file storing editing sites as 'ChrNumber_ChrPos\n'
    :param chr_selected:    chromosome number for which editing sites are to be loaded

    :return: a hash storing editing sites as 'ChrNumber_ChrPos'

    This function loads editing sites from a editing sites file for a selected chromosome number.
    """

    edit_sites_file = open(editsites_file, "r")  # opening editing sites file
    edit_hash = {}  # hash storing editing sites as 'ChrNumber_ChrPos'

    for line in edit_sites_file:  # iterating through editing sites file
        if line.split("_")[0] == chr_selected:  # checking if editing site belongs to selected chromosome
            edit_hash[line.strip()] = 0  # if so, storing editing site to return hash

    edit_sites_file.close()  # closing input editing sites file
    return (edit_hash)  # returning


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

    insertions = list()  # initialize list storing insertion positions (positions w.r.t. reference,
    # position before insertion event, 1-based)
    deletions = list()  # initialize list storing deletion positions (positions w.r.t. reference,
    # position before deletion event, 1-based)
    conf_table = {"AA": 0, "AT": 0, "AG": 0, "AC": 0, "AN": 0, "AI": 0, "AD": 0,
                  "TA": 0, "TT": 0, "TG": 0, "TC": 0, "TN": 0, "TI": 0, "TD": 0,
                  "GA": 0, "GT": 0, "GG": 0, "GC": 0, "GN": 0, "GI": 0, "GD": 0,
                  "CA": 0, "CT": 0, "CG": 0, "CC": 0, "CN": 0, "CI": 0, "CD": 0,  # initialize confusion table
                  "NA": 0, "NT": 0, "NG": 0, "NC": 0, "NN": 0, "NI": 0, "ND": 0,
                  "YA": 0, "YT": 0, "YG": 0, "YC": 0, "YN": 0, "YI": 0, "YD": 0,
                  # added in 2021    # (all counts are 0)
                  "SA": 0, "ST": 0, "SG": 0, "SC": 0, "SN": 0, "SI": 0, "SD": 0,
                  "WA": 0, "WT": 0, "WG": 0, "WC": 0, "WN": 0, "WI": 0, "WD": 0,
                  "KA": 0, "KT": 0, "KG": 0, "KC": 0, "KN": 0, "KI": 0, "KD": 0,
                  "MA": 0, "MT": 0, "MG": 0, "MC": 0, "MN": 0, "MI": 0, "MD": 0,
                  "BA": 0, "BT": 0, "BG": 0, "BC": 0, "BN": 0, "BI": 0, "BD": 0,
                  "DA": 0, "DT": 0, "DG": 0, "DC": 0, "DN": 0, "DI": 0, "DD": 0,
                  "HA": 0, "HT": 0, "HG": 0, "HC": 0, "HN": 0, "HI": 0, "HD": 0,
                  "VA": 0, "VT": 0, "VG": 0, "VC": 0, "VN": 0, "VI": 0, "VD": 0}

    # computing confusion table

    read_pos = 0  # current read position (0-based)
    ref_pos = 0  # current reference position (0-based)
    for c in cigar:  # filling confusion table with the cigar string

        operation = c[0]  # type of operation is encoded by numbers (m, ins, del)
        op_length = c[1]  # length (stretch) of operation
        snp_key_i = chr + "_"  # SNP and A-to-I sites key (initialized, incomplete)

        if operation == 4:
            continue  # skip soft-clipped positions

        if operation == 1:  # CASE: insertion
            ref_nt = ref_seq[ref_pos]  # reference nucleotide
            table_key = ref_nt + "I"  # generate confusion table key
            conf_table[table_key] += op_length  # incrementing insertion counter
            insertions.extend([ref_pos] * op_length)  # storing insertion positions
            read_pos += op_length  # update read position (ref position does not change)

        elif operation == 2:  # CASE: deletion
            ref_nt = ref_seq[ref_pos]  # reference nucleotide
            table_key = ref_nt + "D"  # generate confusion table key
            conf_table[table_key] += op_length  # incrementing deletion counter
            deletions.extend([ref_pos] * op_length)  # storing deletion positions
            ref_pos += op_length  # update reference position (read pos. does not change)

        else:  # CASE: match (real match or mismatch)
            for readp, refp in zip(range(read_pos, read_pos + op_length), range(ref_pos, ref_pos + op_length)):
                if (snp_key_i + str(ref_start + 1 + refp)) in snp_ai_hash:
                    continue  # skipping SNP / A-to-I sites
                if ref_seq[refp] == read_seq[readp]:  # real match
                    table_key = ref_seq[refp] * 2  # generate confusion table key
                    conf_table[table_key] += 1  # incrementing match counter
                else:  # mismatch
                    ref_nt = ref_seq[refp]  # reference nucleotide
                    table_key = ref_nt + read_seq[readp]  # generate confusion table key
                    conf_table[table_key] += 1  # incrementing mismatch counter

            read_pos += op_length  # update read position
            ref_pos += op_length  # update reference position

    return ([insertions, deletions, conf_table])  # returning


def update_region_borders(region, minpos_hash, maxpos_hash, ref_start, ref_end):
    """
    :param region:          genomic region name
    :param minpos_hash:     hash storing region_name -> minimum position covered (0-based, including)
    :param maxpos_hash:     hash storing region_name -> maximum position covered (0-based, excluding)
    :param ref_start:       read starting position
    :param ref_end:         read ending position

    :return:    void

    This function updates the minimum and maximum position of a genomic region covered by reads.
    """

    if minpos_hash[region] > ref_start:
        minpos_hash[region] = ref_start
    if maxpos_hash[region] < ref_end:
        maxpos_hash[region] = ref_end

    return ()


def labeling_ntpair(region_strand):
    """
    :param region_strand:   genomic region strand (':': unspecified, '+': reference reverse, '-': reference forward)
    :return:                a list containing the converted (original) and conversion (observed) nucleotide
    """

    if region_strand == "-":
        converted_nt = "A"  # reference forward region: "A" to be converted (labeled)
        conversion_nt = "G"  # reference forward region: "G" to be observed
    elif region_strand == "+":
        converted_nt = "T"  # reference reverse region: "T" to be converted (labeled)
        conversion_nt = "C"  # reference forward region: "C" to be observed
    else:
        converted_nt = None
        conversion_nt = None

    return ([converted_nt, conversion_nt])


def process_cigar_snps_conv(mode, cigar, converted_nt, conversion_nt,
                            read_id, read_seq, read_chr, read_region, ref_seq, ref_start,
                            region_insinfo, regions_convpos, snps={}):
    """
    :param mode:                mode to be executed; choose 'snp_all' to compute SNPs reporting all kind of conversions,
                                choose 'snp_exclusive' to compute SNPs reporting all but labeling-specific conversions,
                                choose 'conv' to compute conversions only, choose 'both' to compute both SNPs (exclu-
                                sive) and conversions
    :param cigar:               read cigar string to be processed
    :param converted_nt:        nucleotide converted due to metabolic labeling
    :param conversion_nt:       nucleotide observed due to metabolic labeling
    :param read_id:             read ID
    :param read_seq:            read sequence
    :param read_chr:            chromosome number the read is mapped to
    :param read_region:         genomic region name the read is mapped to
    :param ref_seq:             nucleotide sequence the read is mapped to
    :param ref_start:           reference starting position of the read's mapping
    :param region_insinfo:      hash storing read_region -> hash of potential conversion positions (0-based)
    :param regions_convpos:     hash storing read_region -> {position -> {read_ids with insertion at that position}}
    :param snps:                hash storing SNPs (will not be considered for conversions) (default: empty hash)

    :return:    a list containing (in the following order) a list storing read SNP information and a list storing read
                conversion information

    This function extracts read SNP and conversion information from a given cigar string. SNPs are denoted as None if
    the position is not covered by the read, 0 in case no SNP is present, and the SNP nucleotide otherwise (deletions
    are denoted with 'D', insertions with 'I' in additional columns at the end of the table). Conversions are stored as
    None in case that the position is not covered by the read, masked by a deletion or SNP correction, 0 in case no
    conversion is present, and as 1 otherwise.
    """

    read_snps = []  # list storing read SNP information
    read_conv = []  # list storing read conversion information
    read_pos = 0  # current read position (0-based)
    ref_pos = 0  # current reference position (0-based)

    if mode == "snp_all":  # MODE: all SNPs ...................................................................

        for c in cigar:  # filling read SNPs and conversions with the cigar string

            operation = c[0]  # type of operation is encoded by numbers (m, ins, del)
            op_length = c[1]  # length (stretch) of operation

            if operation == 4:
                continue  # skip soft-clipped positions

            if operation == 1:  # --- CASE: insertion ---

                ins_position = ref_start + ref_pos  # insertion position start
                if ins_position in region_insinfo[read_region]:  # storing insertions
                    region_insinfo[read_region][ins_position][read_id] = None
                else:
                    region_insinfo[read_region][ins_position] = {read_id: None}
                read_pos += op_length  # update read position (ref position does not change)

            elif operation == 2:  # --- CASE: deletion ---

                read_snps.extend(["D"] * op_length)  # storing deletions
                ref_pos += op_length  # update reference position (read pos. does not change)

            else:  # --- CASE: match (real match or mismatch) ---

                for readp, refp in zip(range(read_pos, read_pos + op_length), range(ref_pos, ref_pos + op_length)):

                    ref_nt = ref_seq[refp]  # reference sequence nucleotide
                    read_nt = read_seq[readp]  # read nucleotide
                    if ref_nt == read_nt:  # real match
                        read_snps.append(0)  # no SNP recorded
                    else:  # mismatch
                        read_snps.append(read_nt)  # SNP recorded

                read_pos += op_length  # update read position
                ref_pos += op_length  # update reference position

    elif mode == "snp_exclusive":  # MODE: all SNPs but labeling conversions ..........................................

        for c in cigar:  # filling read SNPs and conversions with the cigar string

            operation = c[0]  # type of operation is encoded by numbers (m, ins, del)
            op_length = c[1]  # length (stretch) of operation

            if operation == 4:
                continue  # skip soft-clipped positions

            if operation == 1:  # --- CASE: insertion ---

                ins_position = ref_start + ref_pos  # insertion position start
                if ins_position in region_insinfo[read_region]:  # storing insertions
                    region_insinfo[read_region][ins_position][read_id] = None
                else:
                    region_insinfo[read_region][ins_position] = {read_id: None}
                read_pos += op_length  # update read position (ref position does not change)

            elif operation == 2:  # --- CASE: deletion ---

                read_snps.extend(["D"] * op_length)  # storing deletions
                ref_pos += op_length  # update reference position (read pos. does not change)

            else:  # --- CASE: match (real match or mismatch) ---

                for readp, refp in zip(range(read_pos, read_pos + op_length), range(ref_pos, ref_pos + op_length)):

                    ref_nt = ref_seq[refp]  # reference sequence nucleotide
                    read_nt = read_seq[readp]  # read nucleotide

                    if ref_nt == read_nt:  # real match
                        read_snps.append(0)  # no SNP recorded

                    else:  # mismatch
                        if ref_nt == converted_nt and read_nt == conversion_nt:  # excluding labeling conversions
                            read_snps.append(0)  # no SNP recorded
                        else:
                            read_snps.append(read_nt)  # SNP recorded

                read_pos += op_length  # update read position
                ref_pos += op_length  # update reference position

    elif mode == "conv":  # MODE: conversions ................................................................

        for c in cigar:  # filling read SNPs and conversions with the cigar string

            operation = c[0]  # type of operation is encoded by numbers (m, ins, del)
            op_length = c[1]  # length (stretch) of operation
            snp_key_i = read_chr + "_"  # SNP sites key (initialized, incomplete)

            if operation == 4:
                continue  # skip soft-clipped positions

            if operation == 1:  # --- CASE: insertion ---

                read_pos += op_length  # update read position (ref position does not change)

            elif operation == 2:  # --- CASE: deletion ---

                ref_nts = ref_seq[ref_pos:ref_pos + op_length]  # getting reference nucleotides
                for p, nt in enumerate(ref_nts):  # iterating through reference nucleotides
                    if nt == converted_nt:  # conversion is stored as None (masked by deletion)
                        read_conv.append(None)  # position is stored to position hash
                        regions_convpos[read_region][ref_start + ref_pos + p] = None
                ref_pos += op_length  # update reference position (read pos. does not change)

            else:  # --- CASE: match (real match or mismatch) ---

                for readp, refp in zip(range(read_pos, read_pos + op_length), range(ref_pos, ref_pos + op_length)):

                    ref_nt = ref_seq[refp]  # reference sequence nucleotide
                    read_nt = read_seq[readp]  # read nucleotide

                    if ref_nt == read_nt:  # real match
                        if ref_nt == converted_nt:
                            if (snp_key_i + str(ref_start + 1 + refp)) in snps:  # masked by SNP correction
                                read_conv.append(None)  # (stored as None)
                            else:  # not converted
                                read_conv.append(0)  # (stored as 0)
                            regions_convpos[read_region][ref_start + refp] = None

                    else:  # mismatch
                        if ref_nt == converted_nt:
                            if (snp_key_i + str(ref_start + 1 + refp)) in snps:  # masked by SNP correction
                                read_conv.append(None)  # (stored as None)
                            elif read_nt == conversion_nt:  # is converted
                                read_conv.append(1)  # (stored as 1)
                            else:  # another SNP
                                read_conv.append(0)  # (stored as 0)
                            regions_convpos[read_region][ref_start + refp] = None

                read_pos += op_length  # update read position
                ref_pos += op_length  # update reference position

    elif mode == "both":  # MODE: both SNPs and conversions ..................................................

        for c in cigar:  # filling read SNPs and conversions with the cigar string

            operation = c[0]  # type of operation is encoded by numbers (m, ins, del)
            op_length = c[1]  # length (stretch) of operation
            snp_key_i = read_chr + "_"  # SNP sites key (initialized, incomplete)

            if operation == 4:
                continue  # skip soft-clipped positions

            if operation == 1:  # --- CASE: insertion ---

                ins_position = ref_start + ref_pos  # insertion position start
                if ins_position in region_insinfo[read_region]:  # storing insertions
                    region_insinfo[read_region][ins_position][read_id] = None
                else:
                    region_insinfo[read_region][ins_position] = {read_id: None}
                read_pos += op_length  # update read position (ref position does not change)

            elif operation == 2:  # --- CASE: deletion ---

                ref_nts = ref_seq[ref_pos:ref_pos + op_length]  # getting reference nucleotides
                for p, nt in enumerate(ref_nts):  # iterating through reference nucleotides
                    if nt == converted_nt:  # conversion is stored as None (masked by deletion)
                        read_conv.append(None)  # position is stored to position hash
                        regions_convpos[read_region][ref_start + ref_pos + p] = None
                read_snps.extend(["D"] * op_length)  # storing deletions
                ref_pos += op_length  # update reference position (read pos. does not change)

            else:  # --- CASE: match (real match or mismatch) ---

                for readp, refp in zip(range(read_pos, read_pos + op_length), range(ref_pos, ref_pos + op_length)):

                    ref_nt = ref_seq[refp]  # reference sequence nucleotide
                    read_nt = read_seq[readp]  # read nucleotide

                    if ref_nt == read_nt:  # real match
                        if ref_nt == converted_nt:
                            if (snp_key_i + str(ref_start + 1 + refp)) in snps:  # masked by SNP correction
                                read_conv.append(None)  # (stored as None)
                            else:  # not converted
                                read_conv.append(0)  # (stored as 0)
                            regions_convpos[read_region][ref_start + refp] = None
                        read_snps.append(0)  # no SNP recorded

                    else:  # mismatch
                        if ref_nt == converted_nt:
                            if (snp_key_i + str(ref_start + 1 + refp)) in snps:  # masked by SNP correction
                                read_conv.append(None)  # (stored as None)
                            elif read_nt == conversion_nt:  # is converted
                                read_conv.append(1)  # (stored as 1)
                            else:  # another SNP
                                read_conv.append(0)  # (stored as 0)
                            regions_convpos[read_region][ref_start + refp] = None
                        read_snps.append(read_nt)  # SNP recorded

                read_pos += op_length  # update read position
                ref_pos += op_length  # update reference position

    return (
    [read_snps, read_conv])  # RETURNING ........................................................................


def construct_snp_table(reads, read_snpinfo, region_insinfo, minpos, maxpos):
    """
    :param reads:               list of read IDs to be processed
    :param read_snpinfo:        hash storing read_ID -> [starting position (0-based, incl), SNPs per position]
    :param region_insinfo:      hash storing reference_position -> {read_ids with ins at that position}
    :param minpos:              minimum reference position covered by a read
    :param maxpos:              maximum reference position covered by a read

    :return:    a detailed SNP table; nested list, first index corresponding to rows and second to columns:
                    row 1: genomic region's reference sequence positions (0-based), other rows: reads,
                    column 1: read names, other columns: reference sequence positions,
                    entries: None if the position is not covered by the read, 0 if read has no SNP at that position,
                        SNP nucleotide otherwise

    This function computes a detailed SNP table for all reads of a certain genomic region.
    """

    snptable = []  # initializing SNP table

    region_length = maxpos - minpos  # genomic region total length
    inspos = list(region_insinfo.keys())  # genomic region insertion positions
    inspos.sort()  # sorting insertion positions
    snptable.append([None] + list(range(minpos, maxpos)) + inspos)  # SNP table's first row
    for r in reads:  # iterating through genomic region's reads
        r_start = read_snpinfo[r][0]  # read starting position
        r_entry = [r] + [None] * (r_start - minpos)  # uncovered leading positions are set to None
        r_entry += read_snpinfo[r][1:]  # adding read's SNP information
        r_entry += [None] * (region_length - len(r_entry) + 1)  # uncovered trailing positions are set to None
        for insp in inspos:  # finally adding insertions, too
            if r in region_insinfo[insp]:
                r_entry.append(1)
            else:
                r_entry.append(0)
        snptable.append(r_entry)  # adding read's entry to SNP table

    return (snptable)  # returning SNP table


def construct_conv_table(reads, read_convinfo, region_convpos):
    """
    :param reads:               list of read IDs to be processed
    :param read_convinfo:       hash storing read_ID -> [starting position (0-based, incl), conv per position]
    :param region_convpos:      list of potential conversion positions (0-based)

    :return:    a detailed conversion table; nested list, first index corresponding to rows and second to columns):
                    row 1: potential reference sequence conversion position, other rows: reads,
                    column 1: read names, other columns: potential reference sequence conversion positions,
                    entries: None if the position is not covered, masked by a deletion or SNP correction, 0 if read
                        has no conversion at that position (may have any other SNP, though), 1 otherwise

    This function computes a detailed conversion table for all reads of a certain genomic region.
    """

    convtable = []  # initializing conversion table

    region_convpos.sort()  # sorting potential conversion positions
    conv_length = len(region_convpos)  # number of potential conversion positions
    convtable.append([None] + region_convpos)  # conversion table's first row

    for r in reads:  # iterating through genomic region's reads
        r_start = read_convinfo[r][0]  # read starting positions
        r_entry = [r]  # read ID
        for p in region_convpos:  # uncovered leading positions are set to None
            if p < r_start:
                r_entry.append(None)
            else:
                break
        r_entry += read_convinfo[r][1:]  # adding read's conversion information
        r_entry += [None] * (conv_length - len(r_entry) + 1)  # uncovered trailing positions are set to None
        convtable.append(r_entry)  # adding read's entry to conversion table

    return (convtable)  # returning conversion table


def construct_snp_conv_summary_table(mode, snptable, convtable):
    """
    :param mode:            mode to be executed; choose 'snp' for computing a SNP table summary, choose 'conv' for
                            computing a conversion table summary, choose 'both' for computing both summary tables
    :param snptable:        detailed SNP table; nested list, first index corresponding to rows and second to columns:
                                row 1: genomic region's reference sequence positions (0-based), other rows: reads,
                                column 1: read names, other columns: reference sequence positions,
                                entries: None if the position is not covered by the read, 0 if read has no SNP at that
                                    position, SNP nucleotide otherwise
    :param convtable:       detailed conversion table; nested list, first index corresponding to rows and second to
                            columns):
                                row 1: potential reference sequence conversion position, other rows: reads,
                                column 1: read names, other columns: potential reference sequence conversion positions,
                                entries: None if the position is not covered, masked by a deletion or SNP correction,
                                    0 if read has no conversion at that position (may have any other SNP, though),
                                    1 otherwise

    :return:    a list containing the SNP and conversion summary tables; summary tables are nested lists, with first
                index corresponding to rows and second index to columns; first row stores positions (first entry: None),
                second row stores total counts (first entry: 'total_counts'), third row stores event counts (first
                entry: 'event_counts')

    This function computes summary tables of detailed SNP and conversion tables.
    """

    snptable_summary = []
    convtable_summary = []

    if mode in ("snp", "both"):  # SNP summary

        total_counts = [0 for i in snptable[0][1:]]  # initializing list storing total counts
        event_counts = [0 for i in snptable[0][1:]]  # initializing list storing event counts
        for read in snptable[1:]:  # iterating through reads
            for p, entry in enumerate(read[1:]):  # iterating through read positions
                if entry not in (None, "I"):  # counting up total counts
                    total_counts[p] += 1
                    if entry != 0:  # counting up event counts
                        event_counts[p] += 1
        snptable_summary = [snptable[0], ["total_counts"] + total_counts, ["event_counts"] + event_counts]

    if mode in ("conv", "both"):  # conversion summary

        total_counts = [0 for i in convtable[0][1:]]  # initializing list storing total counts
        event_counts = [0 for i in convtable[0][1:]]  # initializing list storing event counts
        for read in convtable[1:]:  # iterating through reads
            for p, entry in enumerate(read[1:]):  # iterating through read positions
                if entry is not None:  # counting up total counts
                    total_counts[p] += 1
                    if entry == 1:  # counting up event counts
                        event_counts[p] += 1
        convtable_summary = [convtable[0], ["total_counts"] + total_counts, ["event_counts"] + event_counts]

    return ([snptable_summary, convtable_summary])


def write_read_table(r_table, outfile):
    """
    :param r_table:     a hash assigning each genomic region's name a list containing a SNP and a conversion table
                        (implemented as nested lists, the first index corresponding to rows and the second to columns):
                            SNP table: row 1: genomic region's reference sequence positions (0-based), other rows:
                                reads, column 1: read names, other columns: reference sequence positions, entries: None
                                if the position is not covered by the read, 0 if read has no SNP at that position,
                                SNP nucleotide otherwise
                            conversion table: row 1: potential reference sequence conversion position, other rows:
                                reads, column 1: read names, other columns: potential reference sequence conversion
                                positions, entries: None if the position is not covered, masked by a deletion, a SNP or
                                SNP correction, 0 if read has no conversion at that position, 1 otherwise
    :param outfile:     output file to store read tables to

    :return: void

    ....................................................................................................................

    This function writes the contents of a read table hash to a file. The file is formatted as described in the
    following:

    #genomic_region_name\n
    >SNP_TABLE\n
    tab-separated row 1 entries\n
    tab-separated row 2 entries\n
    ...
    >CONV_TABLE\n
    tab-separated row 1 entries\n
    tab-separated row 2 entries\n
    ...
    """

    out_file = open(outfile, "w")  # opening output file
    for gr in r_table:  # iterating through genomic regions
        gr_snps = r_table[gr][0]  # getting genomic region's SNP table
        gr_conv = r_table[gr][1]  # getting genomic region's conversion table
        out_file.write("#" + gr + "\n")  # writing genomic region's name
        out_file.write(">SNP_TABLE\n")  # indicating SNP table listing in the following
        for row in gr_snps:  # iterating through SNP table rows
            for entry in row[:-1]:  # iterating through row entries
                out_file.write(str(entry) + "\t")  # writing row entries (tab-separated)
            out_file.write(str(row[-1]) + "\n")
        out_file.write(">CONV_TABLE\n")  # indicating conversion table listing in the following
        for row in gr_conv:  # iterating through conversion table rows
            for entry in row[:-1]:  # iterating through row entries
                out_file.write(str(entry) + "\t")  # writing row entries (tab-separated)
            out_file.write(str(row[-1]) + "\n")
        out_file.write("\n")

    out_file.close()  # closing output file
    return ()  # returning


def write_edit_sites(rtable_file, editing_cutoff, outfile, chr_numbers=chr_numbers_hs):
    """
    :param rtable_file:     read table file
    :param editing_cutoff:  maximum editing rate (SNP rate) of a nucleotide position to be tolerated
    :param outfile:         output file to write editing sites to

    :param chr_numbers:     hash storing valid chromosome numbers (default: hs chromosome numbers)

    :return: void

    This function detects editing sites, computing position-wise SNP rates from a read table file and storing all
    positions with a SNP rate above a certain threshold. Editing positions are written to an output file as
    'ChrNumber_ChrPos\n'
    """

    rtable = load_read_table(rtable_file)  # loading read table
    out_file = open(outfile, "w")  # opening output file

    for gr in rtable:  # iterating through genomic regions

        gr_name_fields = gr.split("_")  # splitting genomic region name into single parts
        for field in gr_name_fields:  # getting genomic region's corresponding chromosome number
            if field.startswith("chr"):
                gr_chr = field[3:]
        if gr_chr not in chr_numbers:
            continue  # skipping genomic regions not matching a valid chromosome number

        gr_snp_table = rtable[gr][0]  # genomic region's SNP table

        for i, pos in enumerate(gr_snp_table[0][1:]):  # iterating through genomic region positions

            total = gr_snp_table[1][i + 1]  # getting total counts
            if total == 0:
                continue  # (skipping uncovered positions)
            events = gr_snp_table[2][i + 1]  # getting event (SNP) counts
            if (events / total) > editing_cutoff:  # computing SNP rate
                out_file.write(gr_chr + "_" + str(pos + 1) + "\n")  # storing positions with SNP rate above threshold

    out_file.close()  # closing output file
    return ()  # returning


def load_read_table(infile):
    """
    :param infile:  file storing read tables

    :return: read table; a hash assigning each genomic region's name a list containing a SNP and a conversion table
                        (implemented as nested lists, the first index corresponding to rows and the second to columns):
                            SNP table: row 1: genomic region's reference sequence positions (0-based), other rows:
                                reads, column 1: read names, other columns: reference sequence positions, entries: None
                                if the position is not covered by the read, 0 if read has no SNP at that position,
                                SNP nucleotide otherwise
                            conversion table: row 1: potential reference sequence conversion position, other rows:
                                reads, column 1: read names, other columns: potential reference sequence conversion
                                positions, entries: None if the position is not covered, masked by a deletion, a SNP or
                                SNP correction, 0 if read has no conversion at that position, 1 otherwise

    ....................................................................................................................

    This function loads a read table from a file. The file format should be as described in the following:

    #genomic_region_name\n
    >SNP_TABLE\n
    tab-separated row 1 entries\n
    tab-separated row 2 entries\n
    ...
    >CONV_TABLE
    tab-separated row 1 entries\n
    tab-separated row 2 entries\n
    ...
    """

    r_table = {}  # initializing read table (empty hash)
    in_file = open(infile, "r")  # opening input file

    gr = next(in_file).strip().split("#")[1]  # initialize current genomic region (first genomic region in file)
    r_table[gr] = []  # initialize current genomic region's entry (empty list)
    table = []  # initialize current table (empty list)
    for line in in_file:  # iterating through file
        if line.startswith("#"):  # next genomic region encountered
            r_table[gr].append(table)  # storing current conversion table to current genomic region's entry
            table = []  # re-setting table
            gr = line.strip().split("#")[1]  # setting next genomic region
            r_table[gr] = []  # initializing next genomic region's entry
        elif line.startswith(">SNP"):  # SNP table encountered
            pass  # nothing to do, everything was prepared in the step before
        elif line.startswith(">CONV"):  # conversion table encountered
            r_table[gr].append(table)  # storing current SNP table to current genomic region's entry
            table = []  # re-setting table
        else:  # table row
            entries = line.strip().split("\t")  # retrieving single row entries
            # re-transforming entries into integers, Nones and strings
            entries = [int(i) if str_is_int(i) else None if i == 'None' else i for i in entries]
            table.append(entries)  # storing row entries to current table
    r_table[gr].append(table)  # finally, adding last conversion table to last genomic region

    return r_table  # returning read table


def str_is_int(s):
    """
    :param s:   string to be checked whether it represents an integer or not
    :return:    True if the input string represents an integer, False otherwise
    """

    try:
        int(s)
        return (True)
    except ValueError:
        return (False)

# Wrapper --------------------------------------------------------------------------------------------------------------

def create_editfile(inbam, ref, bed, snp, readtable, editsites, readfilter):
    print(f"Creating Read Table for {os.path.basename(inbam)}...", file=sys.stderr)
    read_table_file = read_table('snp_all_summary', bamfile=inbam, reffile=ref, bedfile=bed, snpfile=snp, coverage=20,
                                 refname_to_chr=r_to_c_hg38, readfilter=readfilter)
    write_read_table(r_table=read_table_file, outfile=readtable)
    write_edit_sites(readtable, 0.05, editsites)

    return None


def create_complete_editingfile(snpfile, editsitefile, chromosome):
    """
    This function generates a complete editing file, comprising snps, a-to-i editing sites, and the customly
    defined editing sites. This dictionary will be used later for filtering the reads.
    :param editfile:
    :param snpfile:
    :param aifile:
    :param refname_to_chr:
    :return:
    """
    chromosome = str(chromosome)
    snps_ai = dict()
    if editsitefile is None:
        snps_ai.update(load_snps(snpfile, chromosome))
    if editsitefile is not None:
        snps_ai.update(load_edit_sites(editsitefile, chromosome))
    return snps_ai

alleditsites = dict()
def callresult(editsite):
    global alleditsites
    alleditsites.update(editsite)

def apply_async_with_callback(snpfile, editsites, refname_to_chr, editfile, threads):
    pool = mp.Pool(int(threads))
    for chr in refname_to_chr.values():
        pool.apply_async(create_complete_editingfile, args=(snpfile, editsites, chr), callback=callresult)
    pool.close()
    pool.join()
    with open(editfile, 'wb') as outputfile:
        pickle.dump(alleditsites, outputfile)


def PreProcessReads(
        control,
        output,
        reffile,
        bedfile,
        snpfile,
        threads=1,
        refname_to_chr=r_to_c_hg38,
        readfilter='singleend'
):

    bam = f"{output}filtering/{control}_mapped_filtered.bam"
    readtable = f"{output}readpreprocess/{control}_readtable.txt"
    editsites = f"{output}readpreprocess/{control}_editfile.txt"
    editfile = f"{output}readpreprocess/{control}_alleditsites.pkl"

    print("Detecting Editing Sites...", file=sys.stderr)
    create_editfile(
        inbam=bam,
        ref=reffile,
        bed=bedfile,
        snp=snpfile,
        readtable=readtable,
        editsites=editsites,
        readfilter=readfilter
    )

    print(f"Creating Editing Site File with {int(threads)} Cores...", file=sys.stderr)
    apply_async_with_callback(snpfile, editsites, refname_to_chr, editfile, threads)

    return None
