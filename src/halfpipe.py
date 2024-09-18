import argparse
from createsummary import CreateSummaryFiles
import csv
from filtering import FilterReads, plot_filteredreads
from mapping import MapReads, plot_mappedreads
import os
import pandas as pd
from parameterfit import ParameterFit
from ratioestimation import EstimateRatios
from readpreprocessing import PreProcessReads
import sys
import yaml

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
    filename, fileextension = os.path.splitext(input)
    if fileextension == '.tsv':
        filterstrategy = config["params"]["filterstrategy"]
        if filterstrategy == "pairedend" or filterstrategy == "pseudosingleend":
            tsv_file = open(input)
            read_tsv = csv.reader(tsv_file, delimiter='\t')
            samples = []
            for pair in read_tsv:
                if len(pair) != 4:
                    raise ValueError(
                        'Input tsv-file is not of format read1.fastg(.gz)\tread2.fastg(.gz)\ttimepoint\toutname.')
                samples.append(pair)
            return samples
        else:
            tsv_file = open(input)
            read_tsv = csv.reader(tsv_file, delimiter='\t')
            samples = []
            for pair in read_tsv:
                if len(pair) != 3:
                    raise ValueError('Input tsv-file is not of format sample.fastg(.gz)\ttimepoint\toutname.')
                samples.append(pair)
            return samples
    else:
        raise ValueError('Input file must be tsv.')

def summary_input(input):
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
    :return:
    """

    # Mapping

    if not os.path.exists(f"{config['output']}/mapping"):
        os.makedirs(f"{config['output']}/mapping")

    info = []
    for sample in samples:
        MapReads(
            sample=sample,
            output=config['output'],
            refgenome=config['input']['refgenome'],
            threads=config['params']['cores'],
            gpus=config['params']['gpus'],
            readfilter=config['params']['filterstrategy']
        )
        info.append(sample[2])

    if config['params']['filterstrategy'] == 'pairedend' or config['params']['filterstrategy'] == 'pseudosingleend':
        bamfiles = [f"{config['output']}/mapping/{sample[3]}_mapped.bam" for sample in samples]
    else:
        bamfiles = [f"{config['output']}/mapping/{sample[2]}_mapped.bam" for sample in samples]
    plot_mappedreads(bamfiles=bamfiles, labels=info, output=config['output'])

    # Filtering

    if not os.path.exists(f"{config['output']}/filtering"):
        os.makedirs(f"{config['output']}/filtering")
    output_filtering = f"{config['output']}/filtering/"



    for bam in bamfiles:
        FilterReads(inbam=bam, outbam=output_filtering, readfilter=config['params']['filterstrategy'])

    if config['params']['filterstrategy'] == 'pairedend' or config['params']['filterstrategy'] == 'pseudosingleend':
        bamfiles = [f"{config['output']}/filtering/{sample[3]}_mapped_filtered.bam" for sample in samples]
    else:
        bamfiles = [f"{config['output']}/filtering/{sample[2]}_mapped_filtered.bam" for sample in samples]

    plot_filteredreads(bamfiles=bamfiles, labels=info, output=config['output'])

    return None

def ReadPreProcess():

    if not os.path.exists(f"{config['output']}/readpreprocess"):
        os.makedirs(f"{config['output']}/readpreprocess")

    PreProcessReads(
        control=control,
        output=config['output'],
        reffile=config['input']['refgenome'],
        bedfile=config['input']['bed'],
        snpfile=config['input']['snps'],
        threads=config['params']['cores'],
        readfilter=config['params']['filterstrategy']
    )

    return None

def RatioEstimation():

    if not os.path.exists(f"{config['output']}/ratioestimation"):
        os.makedirs(f"{config['output']}/ratioestimation")

    EstimateRatios(
        samples=samples,
        control=control,
        config=config
    )

    return None

def SummaryFiles(summarysamples):

    if not os.path.exists(f"{config['output']}/summaryfiles"):
        os.makedirs(f"{config['output']}/summaryfiles")

    CreateSummaryFiles(input=summarysamples, config=config)

    return None

def FitParameters(timepoints):

    if not os.path.exists(f"{config['output']}/parameterfit"):
        os.makedirs(f"{config['output']}/parameterfit")
    ParameterFit(f"{config['output']}/summaryfiles/{config['input']['summaryname']}.tsv",
                 f"{config['output']}/parameterfit/testfit.RData", timepoints, config['params']['cores'],
                 config['params']['model'])

# Reading sample info --------------------------------------------------------------------------------------------------

samples = read_input(config['input']['samples'])
if config['params']['filterstrategy'] == 'pairedend' or config['params']['filterstrategy'] == 'pseudosingleend':
    for pair in samples:
        if str(pair[2]) == '0' or pair[2] == 'control':
            control = pair[3]
else:
    for pair in samples:
        if str(pair[1]) == '0' or pair[1] == 'control':
            control = pair[2]


# Running the pipeline -------------------------------------------------------------------------------------------------

def pipe():


    if not os.path.isdir(config['output']):
        os.makedirs(config['output'])

    if args.command == 'pipe':

        if not os.path.isdir(config['output']):
            os.makedirs(config['output'])

        summarysamples = pd.DataFrame(summary_input(config['input']['summmaryfile']))

        if config['params']['model'] == "onecompartment":
            timepoints = summarysamples.iloc[:, 1]
        elif config['params']['model'] == "twocompartment":
            timepoints = summarysamples.iloc[:, 2]
        else:
            raise ValueError("Model must be either 'onecompartment' or 'twocompartment'")

        print("##### PRE-PROCESSING #####\n", file=sys.stdout)
        MapAndFilter()
        ReadPreProcess()

        print("##### MODELING #####\n", file=sys.stdout)
        RatioEstimation()
        SummaryFiles(summarysamples)

        print("##### POST-PROCESSING #####\n", file=sys.stdout)
        FitParameters(timepoints)

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
        summarysamples = pd.DataFrame(summary_input(config['input']['summmaryfile']))
        SummaryFiles(summarysamples)

    if args.command == "fitparameters":
        print("##### RNA HALF-LIFE ESTIMATION #####\n", file=sys.stdout)
        summarysamples = pd.DataFrame(summary_input(config['input']['summmaryfile']))
        if config['params']['model'] == "onecompartment":
            timepoints = summarysamples.iloc[:, 1]
        elif config['params']['model'] == "twocompartment":
            timepoints = summarysamples.iloc[:, 2]
        else:
            raise ValueError("Model must be either 'onecompartment' or 'twocompartment'")
        FitParameters(timepoints)

    return None


if __name__ == '__main__':
    pipe()