#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --output=/lustre/smriti.a/samhita/out/trimmomatic.out # Standard output and error log
#SBATCH --cpus-per-task=10 # Run a task on 10 cpus
#SBATCH --mem=32GB # Use 32GB of memory.

######## TRIMMO0MATIC ######
java -jar /lustre/smriti.a/samhita/tools_softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
/lustre/smriti.a/samhita/samples/IGIB11301413_R1.fastq.gz /lustre/smriti.a/samhita/samples/IGIB11301413_R2.fastq.gz \
/lustre/smriti.a/samhita/results/trimmomatic/IGIB11301413_R1_paired.fq /lustre/smriti.a/samhita/results/trimmomatic/IGIB11301413_R1_unpaired.fq \
/lustre/smriti.a/samhita/results/trimmomatic/IGIB11301413_R2_paired.fq /lustre/smriti.a/samhita/results/trimmomatic/IGIB11301413_R2_unpaired.fq \
ILLUMINACLIP:/lustre/smriti.a/samhita/tools_softwares/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
