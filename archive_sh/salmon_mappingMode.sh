#!/bin/bash
#SBATCH --job-name=salmon_mapping_1413
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --output=/lustre/smriti.a/samhita/out/salmon_mapping_1413.out # Standard output and error log
#SBATCH --cpus-per-task=20 # Run a task on 10 cpus
#SBATCH --mem=32GB # Use 32GB of memory.    

module load gcc-8.4

salmon quant -i /lustre/smriti.a/samhita/reference/salmon_transcripts_index.fa.gz \
-l A \
-1 /lustre/smriti.a/samhita/results/trimmomatic/IGIB11301413_R1_paired.fq \
-2 /lustre/smriti.a/samhita/results/trimmomatic/IGIB11301413_R2_paired.fq \
--validateMappings \
-o /lustre/smriti.a/samhita/results/salmon/salmon_quant_1413_mappingMode