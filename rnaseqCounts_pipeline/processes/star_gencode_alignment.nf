#!/usr/bin/env nextflow

// Process 3b : ALIGNING THE GENOME USING START ALIGNER
process STAR_GENCODE {
    tag "$sample_id"
    //publishDir "${params.outDir}/${sample_id}/star_gencode", mode: 'move'

    publishDir "${params.outDir}/${sample_id}/star_gencode", pattern: "*_Aligned.sortedByCoord.out.bam", mode: 'move'
    publishDir "${params.outDir}/${sample_id}/star_gencode", pattern: "*_Unmapped.out.mate1", mode: 'move'
    publishDir "${params.outDir}/${sample_id}/star_gencode", pattern: "*_Unmapped.out.mate2", mode: 'move'
    publishDir "${params.outDir}/${sample_id}/star_gencode", pattern: "*__starGENOME/*", mode: 'move'
    publishDir "${params.outDir}/${sample_id}/star_gencode", pattern: "*__starpass1/*", mode: 'move'
    publishDir "${params.outDir}/${sample_id}/star_gencode", pattern: "*_Aligned.toTranscriptome.out.bam", mode: 'move'
    publishDir "${params.outDir}/${sample_id}/star_gencode", pattern: "*_Log.final.out", mode: 'move'
    //publishDir "${params.outDir}/${sample_id}/star_gencode", pattern "*_Aligned.toTranscriptome.out.bam", mode: 'move', when: { !workflow.onComplete }
    //publishDir "${params.outDir}/${sample_id}/star_gencode", pattern "*_Log.final.out", mode: 'move', when: { !workflow.onComplete }*/

    input: 
    tuple val(sample_id), path(trimmedReads)
    path(star_index_gencode)

    output: 
    tuple path("*_SJ.out.tab"),
          path("*_Aligned.sortedByCoord.out.bam"),
          path("*_Unmapped.out.mate1"),
          path("*_Unmapped.out.mate2"), emit: starOutput
    //path "${sample_id}_*", mode: 'move'
    path("*_Aligned.toTranscriptome.out.bam"), emit: transcriptomeBAM
    val sample_id, emit: id
    path("*_Log.final.out"), emit: starLog

    script: 
    """
    module load STAR-2.7.8a

    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${star_index_gencode} \
        --readFilesIn ${trimmedReads[0]} ${trimmedReads[1]} \
        --outFileNamePrefix ./${sample_id}_${params.mmNum}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped Fastx \
        --twopassMode Basic \
        --quantMode TranscriptomeSAM \
        --readFilesCommand zcat \
        --outFilterMultimapNmax ${params.mmNum}
    """
}