#!/usr/bin/env nextflow

// Process : RUNNING FEATURE COUNTS
process FEATURE_COUNTS {
    tag "$sample_id"
    label "featureCounts"
    publishDir "${params.outDir}/${sample_id}/featureCounts", mode: 'move'

    input: 
    val(sample_id)
    path(aligned_BAM)
    path(gtf)

    output: 
    path "${sample_id}_featureCounts_herv.*", emit: featureCounts

    script: 
    """
    featureCounts -p -O -T ${task.cpus} \
        -a ${gtf} \
        -o ${sample_id}_featureCounts_herv.txt \
        ${aligned_BAM}
    """
}