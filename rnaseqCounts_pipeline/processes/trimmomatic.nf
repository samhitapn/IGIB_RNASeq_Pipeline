#!/usr/bin/env nextflow

// Process 2 : TRIM THE ADAPTER SEQUNENCES AND LOW QUALITY READS
process TRIM {
    tag "$sample_id"
    publishDir "${params.outDir}/${sample_id}/trimmomatic", mode: 'symlink'

    input: 
    tuple val(sample_id), path(reads)

    output: 
    //val (sample_id), emit: id
    path "${sample_id}_R1_unpaired.fq.gz", emit: R1_up
    path "${sample_id}_R2_unpaired.fq.gz", emit: R2_up
    tuple val(sample_id), path("*_R1_paired.fq.gz"), path("*_R2_paired.fq.gz"), emit: trimOutput
    
    script: 
    """
   java -jar /lustre/smriti.a/samhita/tools_softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
              -threads ${task.cpus} \
              ${reads.join(' ')} \
              ${sample_id}_R1_paired.fq.gz ${sample_id}_R1_unpaired.fq.gz \
              ${sample_id}_R2_paired.fq.gz ${sample_id}_R2_unpaired.fq.gz \
              ILLUMINACLIP:${params.adapters}:2:20:5 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
              
    """
}