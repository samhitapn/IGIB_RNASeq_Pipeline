#!/usr/bin/env nextflow

// Process 3 : RUNNING SALMON IN MAPPING MODE SIRECTLY USING THE TRIMMED FASTQ FILES
process SALMON_MAP_LNCIPEDIA {
    tag "$sample_id"
    publishDir "${params.outDir}/${sample_id}/salmon_map_lncipedia", mode: 'symlink'

    input: 
    tuple val(sample_id), path(trimmedReads)
    path salmon_index_lncipedia

    output: 
    path "${sample_id}_salmon_quant_mappingMode_lncipedia/*"
    path "${sample_id}_salmon_quant_mappingMode_lncipedia/logs/salmon_quant.log", emit salmonMapLog

    script: 
    """
    module load gcc-8.4

    salmon quant -i ${salmon_index_lncipedia} \
                 -l A \
                 -1 ${trimmedReads[0]} \
                 -2 ${trimmedReads[1]} \
                 --validateMappings \
                 --gcBias \
                 --seqBias \
                 -o ${sample_id}_salmon_quant_mappingMode_lncipedia
    """
}