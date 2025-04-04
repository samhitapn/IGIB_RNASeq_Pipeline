#!/bin/bash
#SBATCH --job-name=star_genomeIndex
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --output=/lustre/smriti.a/samhita/out/star_genomeIndex.out # Standard output and error log
#SBATCH --cpus-per-task=40 # Run a task on 10 cpus
#SBATCH --mem=32GB # Use 32GB of memory.    

module load STAR-2.7.8a

STAR \
--runThreadN 40 \
--runMode genomeGenerate \
--genomeDir  /lustre/smriti.a/samhita/reference/STAR_genomeIndices \
--genomeFastaFiles /lustre/smriti.a/samhita/reference/GRCh38.p14.genome.fa \
--sjdbGTFfile /lustre/smriti.a/samhita/reference/gencode.v47.annotation.gtf \
--sjdbOverhang 149

gffread -w ../../reference/salmon_index/salmon_alignmentMode_genome.fa -g /lustre/smriti.a/samhita/reference/GRCh38.p14.genome.fa /lustre/smriti.a/samhita/reference/gencode.v47.annotation.gtf