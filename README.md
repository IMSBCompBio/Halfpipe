# Halfpipe

Halfpipe is designed for the estimation of RNA half-lives from time course data of RNA labeling protocols such as [SLAM-seq](https://www.nature.com/articles/nmeth.4435) or [TimeLapse-seq](https://www.nature.com/articles/nmeth.4582). The pipeline builds on the modeling framework presented in [Müller *et al.* (2024)](https://doi.org/10.1371/journal.pcbi.1012059). To this end, it automates essential pre-processing, model fitting and post-processing steps. Please read our latest [pre-print](https://www.biorxiv.org/content/10.1101/2024.09.19.613510v1) for more information about the underlying methodology of Halfpipe. 



<p align="center">
    <img src="docs/overview.png" width="950"/>
</p>

## Hardware/Software Requirements

Required:

- 64 Bit Linux system 
- x86-64 compatible CPU(s)
- Workload Manager (e.g., [SLURM](https://slurm.schedmd.com/documentation.html))
- [Conda](https://docs.conda.io/projects/conda/en/stable/) Installation

Optional:

- NVIDIA CUDA GPU(s)

## Setup

Simply download this repository to a location of your choice. We recommend using a high performance machine for computation. A sample script `slurm_example.sh` that you can customize to run Halfpipe with the widely used [SLURM workload manager](https://slurm.schedmd.com/documentation.html) is located in the `docs/` folder. 
Next, set up a [conda](https://docs.conda.io/projects/conda/en/stable/) environment using the `halfpipe_condaenv.yml` file in the `config/` folder:

```console
user@foo:~$ conda env create --file=config/halfpipe_condaenv.yml  
```

This yaml-file contains all modules that are necessary to execute Halfpipe. By default, the environment is called `Halfpipe`, but of course you can call it whatever you want. By default, you can activate the environment with:

```console
user@foo:~$ conda activate Halfpipe  
```

In general Halfpipe is invoked with the `halfpipe.py` script:

```console
user@foo:~$ python3 src/halfpipe.py subcommand  
```

A detailed description of the available sub-commands can be found below. But first we want to emphasize the use of the `config.yml` file. This file contains all the parameters that the user needs to specify, and here is a snippet to show 
what this file looks like:

```yaml
---

input:
  samples: path/to/inputsamplesheet.tsv
  bed: path/to/3UTR.bed
  refgenome: path/to/genome.fa
  snps: path/to/SNPfile.vcf
  summaryname: yoursamplename

output: path/to/outputdirectory/ # output directory

...
```

In general, we advise to only work with absolute
instead of relative paths! 

The parameter `samples` defines the path to the sample-sheet (TSV format required!) which lists the corresponding fastq(.gz)-files. You can view sample sheets in the `docs` folder, namely `docs/mocksamplesheet_singleend.tsv` and `docs/mocksamplesheet_pairedend.tsv`. Column 1 (and 2, in the case of paired-end data) defines the path to the fastq(.gz) samples, while the last two columns define the measurement time in minutes and a descriptive name for the sample for later readability.

The parameter “bed” specifies the path to the BED file containing genomically coordinated, annotated 3'UTRs. Such a file can be downloaded from the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables).

The parameter “Refgenome” defines the path to the fasta file of the reference genome. Please use the [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) genome fasta for GRCh38. 

The parameter `snps` defines the path to the SNP annotation file in VCF format. You can download this file from the [NCBI SNPdb](https://www.ncbi.nlm.nih.gov/snp/) database. 

The parameter `summaryname` can be used to specify any name for the resulting summary file. 

The parameter `output` defines the output path. All subfolders and files generated by halfpipe are saved under this path. 

Further parameters such as the library type and the desired model can be configured in the `config.yml` file.

## Commands

### pipe

The `pipe` subcommand executes the entire pipeline for a 4sU time-series experiment. It maps (using NextGenMap) and filters reads from Fastq raw files, corrects the data for potential SNPs and RNA editing sites, and counts the T>C conversion induced by 4sU labeling. Halfpipe uses these counts to estimate the proportions of newly synthesized reads over time (EM algorithm) before fitting either a one- or two-compartment model of RNA metabolism to these quantities.

Usage: 

```console
user@foo:~$ python3 halfpipe.py pipe  
```



### mapandfilter

The `mapandfilter` subcommand maps and filters reads from raw fastq-files using NextGenMap. The output BAM-files can be found in the `mapped/` and `filtered/` subfolders.

```console
user@foo:~$ python3 halfpipe.py mapandfilter  
```



### readpreprocess

The subcommand `readpreprocess` is used to identify putative editing sites (including SNPs, A-to-I editing sites, etc.) from control samples. This step is crucial for the later calculation of the labeling efficiency with the `ratioestimation` subcommand. If this step is not performed, false-positive T>C conversions would bias the downstream analyses. The output files can be found in the `readpreprocess/` subfolder.


```console
user@foo:~$ python3 halfpipe.py readpreprocess  
```
 
### ratioestimation

The subcommand `ratioestimation` is used to estimate the 4sU incorporation rate and accordingly the ratio of new to total RNA for each measured 3'UTR. These quantities are the targets for model fitting. The output files can be found in the subfolder `ratioestimation/`.

Usage: 

```console
user@foo:~$ python3 halfpipe.py ratioestimation  
```

### createsummary

The subcommand `createsummary` summarizes the individual results of the `ratioestimation`. The resulting output file is required for model fitting and can be found in the subfolder `summaryfiles/`.

Usage: 

```console
user@foo:~$ python3 halfpipe.py createsummary  
```

### fitparameters

The subcommand `fitparameters` applies either the one- or two-compartment model to estimate the RNA half-lives. The output files can be found in the `parameterfit/` subfolder.
 

Usage: 

```console
user@foo:~$ python3 halfpipe.py fitparameters  
```
