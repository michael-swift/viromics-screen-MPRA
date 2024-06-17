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

elif [ $1 = "unlock" ]
    then
        snakemake -s $SNAKEFILE $TARGET -F --rerun-incomplete --unlock --cores 1

elif [ $1 = "touch" ]
    then
        snakemake -s $SNAKEFILE $TARGET -F --rerun-incomplete --unlock --touch --cores 1
    
elif [ $1 = "dry" ]
    then
  # Run snakemake
    echo 'running snakemake'
    snakemake --profile profile/ -n

elif [ $1 = "profile" ]
    then
  # Run snakemake
    echo 'running snakemake'
    snakemake --verbose --profile profile/ --software-deployment-method conda apptainer

elif [ $1 = "sbatch" ]
    # Run snakemake as an SBATCH job, good for long workflows (obviously only works on slurm)
    then
    sbatch \
        --ntasks=1 \
        --cpus-per-task=1 \
        --mem=8000 \
        --mail-user=$EMAIL \
        --time 4-0 \
        -p quake \
        -o $SBATCH_LOGFILE \
        -e $SBATCH_LOGFILE_ERR \
        run_snake.sh profile

else
    echo "wrong option"
fi