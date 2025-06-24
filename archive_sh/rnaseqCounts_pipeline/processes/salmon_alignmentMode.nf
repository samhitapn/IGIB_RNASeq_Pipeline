#!/usr/bin/env nextflow

// Process 3 : RUNNING SALMON IN MAPPING MODE SIRECTLY USING THE TRIMMED FASTQ FILES
process SALMON_ALIGN {
    tag "$sample_id"
    label "salmon"
    publishDir "${params.outDir}/${sample_id}/salmon_align", mode: 'move'

    input: 
    val(sample_id)
    path(transcriptome_aligned_BAM)
    path(transcriptome_ref)

    output: 
    path "${sample_id}_salmon_quant_alignmentMode/*", emit: salmonAlignOutput
    //path "temp/*", emit: salmonAlign_MultiQC
    //path "${sample_id}_salmon_quant_alignmentMode/logs/salmon_quant_align.log", emit: salmonAlignLog

    script: 
    """
    module load gcc-8.4

    salmon quant -t ${transcriptome_ref} \
                 -l A \
                 -a ${transcriptome_aligned_BAM} \
                 -o ${sample_id}_salmon_quant_alignmentMode \
                 -p ${task.cpus}
     """
}

/*   #!/bin/bash
    mkdir -p temp
    for file in "${sample_id}_salmon_quant_alignmentMode"/*; do
        filename=\$(basename "\$file")  # Escaped for Nextflow
        cp -r "\$file" "temp/${sample_id}_align_\${filename}"
    done*/