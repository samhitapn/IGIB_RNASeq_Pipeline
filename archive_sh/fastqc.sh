#!/bin/bash
#SBATCH --job-name=fastqc_DEN_afterTrim
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --output=/lustre/smriti.a/samhita/out/fastqc_DEN_afterTrim.out # Standard output and error log
#SBATCH --cpus-per-task=10 # Run a task on 10 cpus
#SBATCH --mem=32GB # Use 32GB of memory.

####### FASTQC   ##########
module load Fastqc

fastqc -o /lustre/smriti.a/samhita/results_ssh/DEN0001/fastqc/ -t 16 /lustre/smriti.a/samhita/results_ssh/DEN0001/trimmomatic/DEN0001_R1_paired.fq /lustre/smriti.a/samhita/results_ssh/DEN0001/trimmomatic/DEN0001_R2_unpaired.fq
