# Benchmarking pipeline to assess the performance of mixmax

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)


The snakemake workflow simulates multiplexing and demultiplexing experiments, uses mixmax to perform sample identification and assesses its performance.


## Usage

To run the pipeline on a HPC-cluster using slurm, run the command

`sbatch --time=24:00:00 --wrap=\"snakemake --cores 10 --configfile config/config_{seed}.yaml --software-deployment-method conda --rerun-incomplete  -p  --keep-going --profile profiles/slurm`

## Dependencies
The pipeline runs with

* `python=3.12`
* `numpy=1.26``
* `pandas=2.1`

All other dependencies will be installed automatically upon first execution of the snakemake pipeline. Note that this requires an internet connection. If internet connection is only available on the login nodes of the HPC cluster, the installation can be run on the login node with the command `snakemake --configfile config/config_{seed}.yaml --software-deployment-method conda --conda-create-envs-only` .