import argparse
from createsummary import CreateSummaryFiles
import csv
from filtering import FilterReads, plot_filteredreads
from mapping import MapReads, plot_mappedreads
import numpy as np
import os
import pandas as pd
from parameterfit import ParameterFit
from ratioestimation import EstimateRatios
from readpreprocessing import PreProcessReads
import sys
import yaml
from misc import r_to_c_m39, r_to_c_hg38
from multiprocessing import set_start_method

# Creating parser ------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(prog='Halfpipe', description='Command line input.',
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
subparsers = parser.add_subparsers(help='sub-command help page.', dest='command')

# HaLfpipe parser ------------------------------------------------------------------------------------------------------

parser_pipe = subparsers.add_parser('pipe', help='pipe help page.')
parser_mapandfilter = subparsers.add_parser('mapandfilter', help='mapandfilter help page.')
parser_readpreprocess = subparsers.add_parser('readpreprocess', help='readpreprocess help page.')
parser_ratioestimation = subparsers.add_parser('ratioestimation', help='ratioestimation help page.')
parser_createsummary = subparsers.add_parser('createsummary', help='createsummary help page.')
parser_fitparameters = subparsers.add_parser('fitparameters', help='fitparameters help page.')
args = parser.parse_args() # initialize parser


# Loading configurations from config.yaml file -------------------------------------------------------------------------

with open("config/config.yml", "r") as configfile:
    config=yaml.load(configfile, Loader=yaml.SafeLoader)

# helper functions -----------------------------------------------------------------------------------------------------
def read_input(input):
    """
    Extract samples according to the inputsamplesheet.tsv file
    """

    filterstrategy = config["params"]["filterstrategy"]
    model = config["params"]["model"] 
    filename, fileextension = os.path.splitext(input)
    samples = []
    timepoints = []
    controls = []

    if fileextension == '.tsv':
        if filterstrategy == "pairedend" or filterstrategy == "pseudosingleend":
            tsv_file = open(input)
            read_tsv = csv.reader(tsv_file, delimiter='\t')
            if model == "twocompartment":
                expectedcolumns = 5 
                ems = 'Input tsv-file is not of format read1.fastq(.gz)\tread2.fastq(.gz)\ttimepoint\toutname\tcompartment.'
            else: 
                expectedcolumns = 4
                ems = 'Input tsv-file is not of format read1.fastq(.gz)\tread2.fastq(.gz)\ttimepoint\toutname.'
            
            for pair in read_tsv:
                if len(pair) != expectedcolumns:
                    raise ValueError(ems)
                if str(pair[2]) == '0' or pair[2] == 'control':
                    controls.append(pair[3])
                samples.append(pair)
                timepoints.append(pair[2])
            
        else:
            tsv_file = open(input)
            read_tsv = csv.reader(tsv_file, delimiter='\t')
            if model == "twocompartment":
                expectedcolumns = 4 
                ems = 'Input tsv-file is not of format read1.fastq(.gz)\ttimepoint\toutname\tcompartment.'
            else: 
                expectedcolumns = 3
                ems = 'Input tsv-file is not of format read1.fastq(.gz)\ttimepoint\toutname.'
            for pair in read_tsv:
                if len(pair) != expectedcolumns:
                    raise ValueError(ems)
                if str(pair[1]) == '0' or pair[1] == 'control':
                    controls.append(pair[2])
                samples.append(pair)
                timepoints.append(pair[1])
    else:
        raise ValueError('Input file must be tsv.')
    
    if model == "twocompartment":
        samples = [samples[:len(samples)//2], samples[len(samples)//2:]] # subdivide into nuclear and cytosolic samples
    else:
        samples = [samples]

    return samples, timepoints, controls
    


def summary_input(input):
    """
    Processes summaryfile.tsv to match samples for downstream modeling purposes.
    """

    filename, fileextension = os.path.splitext(input)
    if fileextension == '.tsv':
        tsv_file = open(input)
        read_tsv = csv.reader(tsv_file, delimiter='\t')
        samples = []
        model = config['params']['model']
        if model == "onecompartment":
            for pair in read_tsv:
                if len(pair) != 2:
                    raise ValueError(
                        'Input tsv-file is not of format sample\ttimepoint')
                if pair[1] == 0 or pair[1] == '0':
                    continue
                samples.append(pair)
        if model == "twocompartment":
            for pair in read_tsv:
                if len(pair) != 3:
                    raise ValueError(
                        'Input tsv-file is not of format nuclearsamples\tcytosolicsamples')
                if pair[2] == 0 or pair[2] == '0':
                    continue
                samples.append(pair)
        return samples
    else:
        raise ValueError('Input file must be tsv.')

def MapAndFilter():
    """
    Wrapper that calls NGM for Mapping and Specific filtering criteria.
    """

    # Mapping

    if not os.path.exists(f"{config['output']}/mapping"):
        os.makedirs(f"{config['output']}/mapping")

    if config['params']['model'] == "twocompartment":
        files = samples[0] + samples[1]
    else: 
        files = samples[0]        

    info = []
    for file in files:
        MapReads(
            sample=file,
            output=config['output'],
            refgenome=config['input']['refgenome'],
            threads=config['params']['cores'],
            gpus=config['params']['gpus'],
            readfilter=config['params']['filterstrategy']
        )
        info.append(file[2])

    if config['params']['filterstrategy'] == 'pairedend' or config['params']['filterstrategy'] == 'pseudosingleend':
        bamfiles = [f"{config['output']}/mapping/{file[3]}_mapped.bam" for file in files]
    else:
        bamfiles = [f"{config['output']}/mapping/{file[2]}_mapped.bam" for file in files]

    plot_mappedreads(bamfiles=bamfiles, labels=info, output=config['output'])

    # Filtering

    if not os.path.exists(f"{config['output']}/filtering"):
        os.makedirs(f"{config['output']}/filtering")
    output_filtering = f"{config['output']}/filtering/"

    org = config['params']['organism']
    if org == "human":
        refname_to_chr = r_to_c_hg38
    elif org == "mouse":
        refname_to_chr = r_to_c_m39
    else:
        raise ValueError("Organism must be 'human' or 'mouse'")

    for bam in bamfiles:
        FilterReads(inbam=bam, outbam=output_filtering, readfilter=config['params']['filterstrategy'], refname_to_chr=refname_to_chr)

    if config['params']['filterstrategy'] == 'pairedend' or config['params']['filterstrategy'] == 'pseudosingleend':
        bamfiles = [f"{config['output']}/filtering/{file[3]}_mapped_filtered.bam" for file in files]
    else:
        bamfiles = [f"{config['output']}/filtering/{file[2]}_mapped_filtered.bam" for file in files]

    plot_filteredreads(bamfiles=bamfiles, labels=info, output=config['output'])

    return None

def ReadPreProcess():
    """
    Wrapper that computes sequencing errors from the sequencing data.
    """

    if not os.path.exists(f"{config['output']}/readpreprocess"):
        os.makedirs(f"{config['output']}/readpreprocess")
    
    for control in controls:
        
        PreProcessReads(
            control=control,
            output=config['output'],
            reffile=config['input']['refgenome'],
            bedfile=config['input']['bed'],
            snpfile=config['input']['snps'],
            threads=config['params']['cores'],
            readfilter=config['params']['filterstrategy'],
            org = config['params']['organism']
        )

    return None

def RatioEstimation():
    """
    Wrapper that estimates the proportion of newly synthesized transcripts per 3'UTR over time
    """

    if not os.path.exists(f"{config['output']}/ratioestimation"):
        os.makedirs(f"{config['output']}/ratioestimation")

    for counter, sample in enumerate(samples):
        EstimateRatios(
            samples=sample,
            control=controls[counter],
            config=config
        )

    return None

def SummaryFiles():
    """
    Wrapper that summarizes the estimation results.
    """

    if not os.path.exists(f"{config['output']}/summaryfiles"):
        os.makedirs(f"{config['output']}/summaryfiles")

    summarysamples = []
    if config['params']['model'] == "twocompartment":
        files = samples[0] + samples[1]
        index = 3
    else: 
        files = samples[0]
        index = 2

    for file in files:  
        summarysamples.append(file[index]) 

    CreateSummaryFiles(input=summarysamples, config=config)

    return None

def FitParameters():
    """
    Wrapper that calls either a one- or two-compartment model to estimate transcript dynamics.
    """

    if not os.path.exists(f"{config['output']}/parameterfit"):
        os.makedirs(f"{config['output']}/parameterfit")

    ParameterFit(
        f"{config['output']}/summaryfiles/{config['input']['summaryname']}_nocontrol.tsv",
        f"{config['output']}/parameterfit/testfit.RData", np.sort(np.unique(timepoints)[1:]), config['params']['cores'],
        config['params']['model']
    )

# Reading sample info --------------------------------------------------------------------------------------------------

samples, timepoints, controls = read_input(config['input']['samples'])

# Running the pipeline -------------------------------------------------------------------------------------------------

def pipe():
    """
    Function that defines the whole pipeline.
    """

    if not os.path.isdir(config['output']):
        os.makedirs(config['output'])

    if args.command == 'pipe':

        if not os.path.isdir(config['output']):
            os.makedirs(config['output'])

        print("##### PRE-PROCESSING #####\n", file=sys.stdout)
        #MapAndFilter()
        #ReadPreProcess()

        print("##### MODELING #####\n", file=sys.stdout)
        RatioEstimation()
        SummaryFiles()

        print("##### POST-PROCESSING #####\n", file=sys.stdout)
        FitParameters()

    if args.command == 'mapandfilter':
        print("##### MAPPING AND FILTERING #####\n", file=sys.stdout)
        MapAndFilter()

    if args.command == "readpreprocess":
        print("##### PRE-PROCESSING #####\n", file=sys.stdout)
        ReadPreProcess()

    if args.command == "ratioestimation":
        print("##### ESTMATING NEW/TOTAL RNA RATIOS #####\n", file=sys.stdout)
        RatioEstimation()

    if args.command == "createsummary":
        print("##### CREATING SUMMARY FILES #####\n", file=sys.stdout)
        SummaryFiles()

    if args.command == "fitparameters":
        print("##### RNA HALF-LIFE ESTIMATION #####\n", file=sys.stdout)
        FitParameters()

    return None


if __name__ == '__main__':
    set_start_method("spawn", force=True)
    pipe() # run and grab a coffee