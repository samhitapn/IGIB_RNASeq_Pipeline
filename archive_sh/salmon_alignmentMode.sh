#!/bin/bash
#SBATCH --job-name=salmon_alignment_1413
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --output=/lustre/smriti.a/samhita/out/salmon_alignment_1413.out # Standard output and error log
#SBATCH --cpus-per-task=20 # Run a task on 10 cpus
#SBATCH --mem=32GB # Use 32GB of memory.    

module load gcc-8.4

salmon quant -t /lustre/smriti.a/samhita/reference/GRCh38.p14.genome.fa \
-l A \
-a /lustre/smriti.a/samhita/results/STAR_alignment/IGIB11301413Aligned.sortedByCoord.out.bam \
--output /lustre/smriti.a/samhita/results/salmon/salmon_quant_1413