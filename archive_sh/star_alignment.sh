#!/bin/bash
#SBATCH --job-name=star_alignment
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --output=/lustre/smriti.a/samhita/out/star_alignemnt_1413.out # Standard output and error log
#SBATCH --cpus-per-task=20 # Run a task on 10 cpus
#SBATCH --mem=32GB # Use 32GB of memory.    

module load STAR-2.7.8a

STAR \
--runThreadN 20 \
--readFilesCommand zcat \
--genomeDir /lustre/smriti.a/samhita/reference/STAR_genomeIndices \
--readFilesIn /lustre/smriti.a/samhita/results/trimmomatic/IGIB11301413_R1_paired.fq /lustre/smriti.a/samhita/results/trimmomatic/IGIB11301413_R2_paired.fq \
--outFileNamePrefix /lustre/smriti.a/samhita/results/STAR_alignment/IGIB11301413 \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--twopassMode Basic \
--quantMode TranscriptomeSAM