#!/bin/bash

#name shell variables for calling in snakemake and sbatch
DATETIME=$(date "+%Y_%m_%d_%H_%M_%S")
RESTART=0

SNAKEFILE=Snakefile

#Snakemake config
NJOBS=20
WAIT=120

# Initialize Conda. Adjust the path to where your conda.sh is located.
source ~/miniconda3/etc/profile.d/conda.sh


## activate conda environment called snakemake
conda activate snakemake

# for a specific cell / output
TARGET=''
if [ $# -eq 0 ]
  then
    # Dry run snakemake
    snakemake -s $SNAKEFILE $TARGET --use-conda --rerun-incomplete -n --keep-going --rerun-triggers mtime
    
elif [ $1 = "dry" ]
    then
  # Run snakemake
    echo 'running snakemake'
    snakemake --profile profile/ -n

## Use this to run the workflow
elif [ $1 = "profile" ]
    then
  # Run snakemake
    echo 'running snakemake'
    snakemake --verbose --profile profile/ --software-deployment-method conda apptainer

else
    echo "wrong option"
fi
