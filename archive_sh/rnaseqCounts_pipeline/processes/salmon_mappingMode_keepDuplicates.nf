#!/usr/bin/env nextflow

// Process 3 : RUNNING SALMON IN MAPPING MODE SIRECTLY USING THE TRIMMED FASTQ FILES
process SALMON_MAP_KD {
    tag "$sample_id"
    label "salmon"
    publishDir "${params.outDir}/${sample_id}/salmon_map_kd", mode: 'move'

    input: 
    tuple val(sample_id), path(trimmedReads)
    path (salmon_index)

    output: 
    path "*_salmon_quant_mappingMode_kd/*", emit: salmonMapOutput
    //path "temp/*", emit: salmonMap_MultiQC
    //path "*_salmon_quant_mappingMode/logs/salmon_quant.log", emit: salmonMapLog

    script: 
    """
    module load gcc-8.4

    salmon quant -i ${salmon_index} \
                 -l A \
                 -1 ${trimmedReads[0]} \
                 -2 ${trimmedReads[1]} \
                 --validateMappings \
                 --gcBias \
                 --seqBias \
                 -o ${sample_id}_salmon_quant_mappingMode_kd \
                 -p ${task.cpus}
          
    """
}

/* 
    #!/bin/bash
    mkdir -p temp
    for file in "${sample_id}_salmon_quant_mappingMode"/*; do
        filename=\$(basename "\$file")  # Escaped for Nextflow
        cp -r "\$file" "temp/${sample_id}_map_\${filename}"
    done   */