#!/bin/bash -l
#SBATCH --account=maxmustermann
#SBATCH --output=/path/to/output.txt
#SBATCH --error=/path/to/error.txt

# Of course you can add more SBATCH commands depending on your system

# Go to Halfpipe dir
cd /path/to/Halfpipe/

# activate conda environment
conda activate Halfpipe

# Run halfpipe an grab a coffee
python3 /path/to/Halfpipe/src/halfpipe.py pipe
