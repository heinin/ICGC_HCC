#!/bin/bash
#SBATCH --job-name=infer_strandedness
#SBATCH --mail-type=END
#SBATCH --mail-user=hnatri@asu.edu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task 1
#SBATCH --time 96:00:00

export PERL5LIB=/packages/6x/vcftools/0.1.12b/lib/perl5/site_perl

newgrp combinedlab

source activate rseqc_environment

/home/hnatri/Tools/RSeQC-2.6.4/scripts/infer_experiment.py -r /home/hnatri/Tools/gencode.v25.chr_patch_hapl_scaff.annotation.bed -i /mnt/storage/CANCER_DOWNLOADS/PROCESSED/ICGC_HCC/processed_bam/${SAMPLE}_HISAT2_sortedbycoord.bam > /mnt/storage/CANCER_DOWNLOADS/PROCESSED/ICGC_HCC/RSeQC_result/${SAMPLE}HISAT2_RSeQC_inferstrandedness.tsv
