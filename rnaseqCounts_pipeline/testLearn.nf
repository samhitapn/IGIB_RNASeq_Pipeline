#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// CHECK FOR MAIN OUTDIR, IF NOT CREATE
def outDir = file(params.outDir)
if (!outDir.exists()) {
    println "Creating output directory: ${params.outDir}"
    outDir.mkdirs()
}

//CREATING A CHANNEL OF PE FILES
rawReadPairs = Channel.fromFilePairs(params.readsDir)

// INCLUDE THE PROCESSES
include { FASTQC } from './processes/fastqc.nf'
include { FASTQC as FASTQC_POST } from './processes/fastqc.nf'
include { TRIM } from './processes/trimmomatic.nf'
include { SALMON_MAP } from './processes/salmon_mappingMode.nf'
include { STAR_GENCODE } from './processes/star_gencode_alignment.nf'
include { SALMON_ALIGN } from './processes/salmon_alignmentMode.nf'
include { MULTI_QC } from './processes/multiqc.nf'

// ALL WORKFLOW -> WF_ALL
workflow {
    
    FASTQC         (rawReadPairs.map { id, reads -> tuple(id, reads, "preTrim") })  
    
    TRIM           (rawReadPairs)
  
    FASTQC_POST    (TRIM.out.trimOutput.map {id,R1,R2 -> tuple(id, [R1, R2], "postTrim")})

    SALMON_MAP     (TRIM.out.trimOutput.map {id,R1,R2 -> tuple(id, [R1, R2])}, params.salmon_index)

    STAR_GENCODE   (TRIM.out.trimOutput.map {id,R1,R2 -> tuple(id, [R1, R2])}, params.star_index_gencode)

    SALMON_ALIGN   (STAR_GENCODE.out.id, STAR_GENCODE.out.transcriptomeBAM, params.transcriptome_ref)

    /**FASTQC.out.zip.mix(FASTQC_POST.out.zip,
                        FASTQC.out.html, FASTQC_POST.out.html,
                        SALMON_MAP.out.salmonMapOutput,
                       STAR_GENCODE.out.starLog,
                        SALMON_ALIGN.out.salmonAlignOutput
                    ).collect().view() **/
    SALMON_MAP.out.salmonMapOutput.view()
    /**
    MULTI_QC       (FASTQC.out.zip.mix(FASTQC_POST.out.zip,
                        FASTQC.out.html, FASTQC_POST.out.html,
                        SALMON_MAP.out.salmonMapOutput,
                        STAR_GENCODE.out.starLog,
                        SALMON_ALIGN.out.salmonAlignOutput
                    ).collect())**/
}

// PSEUDO ALIGNMENT ONLY WORKFLOW
/**workflow WF_QUANT{
    
    FASTQC         (rawReadPairs,"preTrim")

    TRIM           (rawReadPairs)

    FASTQC         (Channel.of([TRIM.out.id,[TRIM.out.R1_paired, TRIM.out.R2_paired]]),"postTrim")

    SALMON_MAP     (TRIM.out.id, TRIM.out.R1_paired, TRIM.out.R2_paired, params.salmon_index)

    //STAR_GENCODE   (TRIM.out.id, TRIM.out.R1_paired, TRIM.out.R2_paired, params.star_index_gencode)

    //SALMON_ALIGN   (STAR_GENCODE.out.id, STAR_GENCODE.out.transcriptomeBAM, params.transcriptome_ref)

    MULTI_QC       ()
}

// STAR + sALMON ONLY WORKFLOW
workflow WF_MAP_QUANT{
    
    FASTQC         (rawReadPairs,"preTrim")

    TRIM           (rawReadPairs)

    FASTQC         (Channel.of([TRIM.out.id,[TRIM.out.R1_paired, TRIM.out.R2_paired]]),"postTrim")

    //SALMON_MAP     (TRIM.out.id, TRIM.out.R1_paired, TRIM.out.R2_paired, params.salmon_index)

    STAR_GENCODE   (TRIM.out.id, TRIM.out.R1_paired, TRIM.out.R2_paired, params.star_index_gencode)

    SALMON_ALIGN   (STAR_GENCODE.out.id, STAR_GENCODE.out.transcriptomeBAM, params.transcriptome_ref)

    MULTI_QC       ()
}**/