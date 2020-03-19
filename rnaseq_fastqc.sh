#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --mail-type=END
#SBATCH --mail-user=hnatri@asu.edu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 1
#SBATCH --time 96:00:00

export PERL5LIB=/packages/6x/vcftools/0.1.12b/lib/perl5/site_perl

newgrp combinedlab

source activate fastqc_environment

snakemake --snakefile rnaseq_qc.snakefile -j 100 --rerun-incomplete --cluster "sbatch -n 1 --nodes 1 -c 2 -t 96:00:00"