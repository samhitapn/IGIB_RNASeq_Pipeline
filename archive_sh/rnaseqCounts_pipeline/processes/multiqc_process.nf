#!/usr/bin/env nextflow

// Process 1 : Quality check with fastqc
process MULTI_QC {
    //tag "${sample_id}"
    publishDir "${params.outDir}/multiqc", mode: 'move'

    input: 
    path '*'
    //val sample_id

    output: 
    path "${params.runName}.html", emit: multiqc
    path "${params.runName}_data"

    script:
    """
    multiqc . -f -n ${params.runName}.html
    """
}