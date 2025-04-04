#!/bin/bash
#SBATCH --job-name=salmon_generateIndex_lncipedia
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --output=/lustre/smriti.a/samhita/out/salmon_generateIndex_lncipedia.out # Standard output and error log
#SBATCH --cpus-per-task=20 # Run a task on 10 cpus
#SBATCH --mem=32GB # Use 32GB of memory.    

module load gcc-8.4

grep "^>" <(gunzip -c /lustre/smriti.a/samhita/reference/GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > /lustre/smriti.a/samhita/reference/salmon_index_lncipedia/salmon_decoys_lncipedia.txt

sed -i.bak -e 's/>//g' /lustre/smriti.a/samhita/reference/salmon_index_lncipedia/salmon_decoys_lncipedia.txt

#gzip /lustre/smriti.a/samhita/reference/LNCipedia/lncipedia_5_1.fasta /lustre/smriti.a/samhita/reference/LNCipedia/lncipedia_5_1.fa.gz

cat /lustre/smriti.a/samhita/reference/LNCipedia/lncipedia_5_1.fa.gz /lustre/smriti.a/samhita/reference/GRCh38.primary_assembly.genome.fa.gz > /lustre/smriti.a/samhita/reference/salmon_index_lncipedia/salmon_transcripts_index_lncipedia.fa.gz

salmon index -t /lustre/smriti.a/samhita/reference/salmon_index_lncipedia/salmon_transcripts_index_lncipedia.fa.gz -d /lustre/smriti.a/samhita/reference/salmon_index_lncipedia/salmon_decoys_lncipedia.txt -p 12 -i /lustre/smriti.a/samhita/reference/salmon_index_lncipedia/ --gencode