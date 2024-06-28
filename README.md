# viromics-screen-MPRA

- added a snakemake workflow which attempts to recapitulate the data pre-processing and ultimately generate the intermediate counts files used in the notebooks

Snakemake Usage:

Need to install snakemake in a conda environment

Need to run snakemake using the wrapper script run_snake.sh:

  bash run_snake.sh

TODO: 
  - download data manually or,
  - programatically download the data from zenodo using the .json file in data/
     - add that download script as a rule to the beginning of the Snakemake
  - change the analysis pipeline to be more in line with the reported analysis in the paper:
      - e.g. they use bowtie for mapping and we use bwa2
  - exploratory data analysis on the pipeline outputs (reads per sample, reads per viral MPRA construct, technical variability between samples
