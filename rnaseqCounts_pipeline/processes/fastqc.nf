#!/usr/bin/env nextflow

// Process : FastQC of all the processes
process FASTQC {
    tag "${sample_id}_${stage}"
    publishDir "${params.outDir}/${sample_id}/fastqc/${sample_id}_${stage}", mode: 'symlink'

    input: 
    tuple val(sample_id), path(reads), val(stage)
    /**val(sample_id)
    path(reads)
    val(stage)**/

    output: 
    //val (sample_id)
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html
 
    script: 
    """
    fastqc -o . ${reads.join(' ')}
    """
}