#!/bin/bash
#SBATCH --job-name=salmon_generateIndex
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --output=/lustre/smriti.a/samhita/out/salmon_generateIndex.out # Standard output and error log
#SBATCH --cpus-per-task=20 # Run a task on 10 cpus
#SBATCH --mem=32GB # Use 32GB of memory.    

module load gcc-8.4

#grep "^>" <(gunzip -c /lustre/smriti.a/samhita/reference/GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > /lustre/smriti.a/samhita/reference/salmon_index/salmon_decoys.txt

#sed -i.bak -e 's/>//g' /lustre/smriti.a/samhita/reference/salmon_index/salmon_decoys.txt

#cat /lustre/smriti.a/samhita/reference/gencode.v47.transcripts.fa.gz /lustre/smriti.a/samhita/reference/GRCh38.primary_assembly.genome.fa.gz > /lustre/smriti.a/samhita/reference/salmon_index/salmon_transcripts_index.fa.gz

salmon index -t /lustre/smriti.a/samhita/reference/salmon_index/salmon_transcripts_index.fa.gz -d /lustre/smriti.a/samhita/reference/salmon_index/salmon_decoys.txt -p 12 -i /lustre/smriti.a/samhita/reference/salmon_index_keepDuplicates/ --gencode --keepDuplicates