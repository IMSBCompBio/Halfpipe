from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pysam
from readpreprocessing import r_to_c_hg38, r_to_c_m39

def FilterReads(
    inbam,
    outbam,
    readfilter="singleend",
    MQ=2,
    minIdentity=0.95,
    NM=-1,
    refname_to_chr=r_to_c_hg38,
):
    """
    This function filters out reads that do not align to an annotated 3'UTR
    """

    outbam = f"{outbam}/{os.path.splitext(os.path.basename(inbam))[0]}_filtered.bam"
    infile = pysam.AlignmentFile(inbam, "rb")
    outfile = pysam.AlignmentFile(outbam, "wb", template=infile)

    if readfilter == "pairedend" or readfilter == "pseudosingleend":
        for read, mate in read_pair_generator(infile):
            if (
                read.reference_name not in refname_to_chr
                or mate.reference_name not in refname_to_chr
            ):
                continue
            if read.is_unmapped or mate.is_unmapped:
                continue
            if read.mapping_quality < MQ or mate.mapping_quality < MQ:
                continue
            if (
                float(read.get_tag("XI")) < minIdentity
                or float(mate.get_tag("XI")) < minIdentity
            ):
                continue
            if (NM > -1 and int(read.get_tag("NM")) > NM) or (
                NM > -1 and int(mate.get_tag("NM")) > NM
            ):
                continue

            outfile.write(read)
            outfile.write(mate)

    elif readfilter == "singleend":
        for read in infile.fetch():
            if read.reference_name not in refname_to_chr:
                continue
            if read.is_unmapped:
                continue
            if read.mapping_quality < MQ:
                continue
            if float(read.get_tag("XI")) < minIdentity:
                continue
            if NM > -1 and int(read.get_tag("NM")) > NM:
                continue

            outfile.write(read)

    infile.close()
    outfile.close()

    pysam.sort("-o", outbam, outbam)
    pysam.index(outbam)

    return None


def plot_filteredreads(
        bamfiles,
          labels,
            output
            ):
    """
    A simple plotting function to visualize the impact of the read filtering.
    """

    mapped_reads = []
    unmapped_reads = []
    ambiguous_reads = []

    for bam in bamfiles:
        bamfile = pysam.AlignmentFile(bam, "rb")
        mapped = 0
        unmapped = 0
        ambiguous = 0
        for read in bamfile:
            if (
                not read.is_unmapped
                and not read.is_secondary
                and not read.is_supplementary
            ):
                mapped += 1
            if read.is_unmapped:
                unmapped += 1
            if read.is_secondary or read.is_supplementary:
                ambiguous += 1
        bamfile.close()
        mapped_reads.append(mapped)
        unmapped_reads.append(unmapped)
        ambiguous_reads.append(ambiguous)
    fig, ax = plt.subplots()
    ax.bar(
        labels, mapped_reads, 0.35, label="mapped", edgecolor="black", color="#999999"
    )
    ax.bar(
        labels,
        unmapped_reads,
        0.35,
        bottom=mapped_reads,
        label="unmapped",
        edgecolor="black",
        color="#E69F00",
    )
    ax.bar(
        labels,
        ambiguous_reads,
        0.35,
        bottom=np.array(mapped_reads) + np.array(unmapped_reads),
        label="ambiguous",
        edgecolor="black",
        color="#56B4E9",
    )
    ax.legend(loc="upper right")
    ax.set_ylabel("# reads")
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="gray", linestyle="dashed")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_xticks(ticks=np.arange(len(labels)), labels=labels, rotation=45, ha="right")
    fig.tight_layout()
    plt.savefig(f"{output}/filtering/filteringstats_barplot.png", dpi=300)
    fig.show()

    data = [labels, mapped_reads, unmapped_reads, ambiguous_reads]
    data = pd.DataFrame(data).transpose()

    data.columns = ["sample", "mapped_reads", "unmapped_reads", "ambiguous_reads"]
    data.to_csv(f"{output}/filtering/filteringstats.csv", header=True, index=None)

    return None

def read_pair_generator(
    bam, 
    region_string=None
):  # copied from Biostars: https://www.biostars.org/p/306041/
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

    return None