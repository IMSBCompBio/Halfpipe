import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pysam


def MapReads(
    sample,
    output, 
    refgenome, 
    threads=None, 
    gpus=None, 
    readfilter="singleend"
):
    """
    This function specifies the necessary parameters to run NGM for sequence alignment of the corresponding SLAM-seq data.
    """

    if readfilter == "pairedend" or readfilter == "pseudosingleend":
        print("Running NGM in paired-end mode...")
        bamfile = f"{output}/mapping/{sample[3]}_mapped.bam"

        if threads is not None and gpus is None:
            task = f"ngm -1 {sample[0]} -2 {sample[1]} -r {refgenome} -o {bamfile} -5 12 -b --slam-seq 2 --max-polya 4 --strata -l --rg-id readgroup --rg-sm NA:NA:NA -t {threads}"

        if threads is None and gpus is not None:
            if len(gpus) > 1:
                arr = list(range(2))
                sarr = [str(a) for a in arr]
                task = f"ngm -1 {sample[0]} -2 {sample[1]} -r {refgenome} -o {bamfile} -5 12 -b --slam-seq 2 --max-polya 4 --strata -l --rg-id readgroup --rg-sm NA:NA:NA -g {','.join(sarr)}"
            else:
                task = f"ngm -1 {sample[0]} -2 {sample[1]} -r {refgenome} -o {bamfile} -5 12 -b --slam-seq 2 --max-polya 4 --strata -l --rg-id readgroup --rg-sm NA:NA:NA -g 0"

        if threads is not None and gpus is not None:
            if len(gpus) > 1:
                arr = list(range(2))
                sarr = [str(a) for a in arr]
                task = f"ngm -1 {sample[0]} -2 {sample[1]} -r {refgenome} -o {bamfile} -5 12 -b --slam-seq 2 --max-polya 4 --strata -l --rg-id readgroup --rg-sm NA:NA:NA -t {threads} -g {','.join(sarr)}"
            else:
                task = f"ngm -1 {sample[0]} -2 {sample[1]} -r {refgenome} -o {bamfile} -5 12 -b --slam-seq 2 --max-polya 4 --strata -l --rg-id readgroup --rg-sm NA:NA:NA -t {threads} -g 0"

        if threads is None and gpus is None:
            task = f"ngm -1 {sample[0]} -2 {sample[1]} -r {refgenome} -o {bamfile} -5 12 -b --slam-seq 2 --max-polya 4 --strata -l --rg-id readgroup --rg-sm NA:NA:NA"

    elif readfilter == "singleend":
        print("Running NGM in single-end mode...")
        bamfile = f"{output}/mapping/{sample[2]}_mapped.bam"

        if threads is not None and gpus is None:
            task = f"ngm -q {sample[0]} -r {refgenome} -o {bamfile} -5 12 -n 100 -b --slam-seq 2 --max-polya 4 --strata -l --rg-id readgroup --rg-sm NA:NA:NA -t {threads}"

        if threads is None and gpus is not None:
            if len(gpus) > 1:
                arr = list(range(2))
                sarr = [str(a) for a in arr]
                task = f"ngm -q {sample[0]} -r {refgenome} -o {bamfile} -5 12 -n 100 -b --slam-seq 2 --max-polya 4 --strata -l --rg-id readgroup --rg-sm NA:NA:NA -g {','.join(sarr)}"
            else:
                task = f"ngm -q {sample[0]} -r {refgenome} -o {bamfile} -5 12 -n 100 -b --slam-seq 2 --max-polya 4 --strata -l --rg-id readgroup --rg-sm NA:NA:NA -g 0"

        if threads is not None and gpus is not None:
            if len(gpus) > 1:
                arr = list(range(2))
                sarr = [str(a) for a in arr]
                task = f"ngm -q {sample[0]} -r {refgenome} -o {bamfile} -5 12 -n 100 -b --slam-seq 2 --max-polya 4 --strata -l --rg-id readgroup --rg-sm NA:NA:NA -t {threads} -g {','.join(sarr)}"
            else:
                task = f"ngm -q {sample[0]} -r {refgenome} -o {bamfile} -5 12 -n 100 -b --slam-seq 2 --max-polya 4 --strata -l --rg-id readgroup --rg-sm NA:NA:NA -t {threads} -g 0"

        if threads is None and gpus is None:
            task = f"ngm -q {sample[0]} -r {refgenome} -o {bamfile} -5 12 -n 100 -b --slam-seq 2 --max-polya 4 --strata -l --rg-id readgroup --rg-sm NA:NA:NA"

    os.system(f"{task} --no-progress")  # Running NGM

    pysam.sort("-o", bamfile, bamfile)
    pysam.index(bamfile)

    return None


def plot_mappedreads(
    bamfiles, 
    labels, 
    output
):
    
    """
    A simple plotting function to visualized the mapping results.
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
    plt.savefig(f"{output}/mapping/mappingstats_barplot.png", dpi=300)
    plt.close()

    data = [labels, mapped_reads, unmapped_reads, ambiguous_reads]
    data = pd.DataFrame(data).transpose()

    data.columns = ["sample", "mapped_reads", "unmapped_reads", "ambiguous_reads"]
    data.to_csv(f"{output}/mapping/mappingstats.csv", header=True, index=None)

    return None
