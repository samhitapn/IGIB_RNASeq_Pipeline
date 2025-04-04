#!/usr/bin/env nextflow
nextflow.enable.dsl=2
//nextflow.preview.output = true

// CHECK FOR MAIN OUTDIR, IF NOT CREATE
def outDir = file(params.outDir)
if (!outDir.exists()) {
  println "Creating output directory: ${params.outDir}"
  outDir.mkdirs()
}

//CREATING A CHANNEL OF PE FILES
trimmedReadPairs = Channel.fromFilePairs(params.trimmedReadsDir)

// INCLUDE THE PROCESSES
//include { FASTQC } from './processes/fastqc.nf'
//include { FASTQC as FASTQC_POST } from './processes/fastqc.nf'
//include { TRIM } from './processes/trimmomatic.nf'
include { SALMON_MAP } from './processes/salmon_mappingMode.nf'
include { STAR_GENCODE } from './processes/star_gencode_alignment.nf'
include { SALMON_ALIGN } from './processes/salmon_alignmentMode.nf'
include { MULTI_QC } from './processes/multiqc.nf'
include { FILE_CLEANUP } from './processes/fileCleanup_star.nf'

// ALL WORKFLOW -> WF_ALL
workflow {
    
    //FASTQC         (rawReadPairs.map { id, reads -> tuple(id, reads, "preTrim") })  
    
    //TRIM           (rawReadPairs)
  
    //FASTQC_POST    (TRIM.out.trimOutput.map {id,R1,R2 -> tuple(id, [R1, R2], "postTrim")})
    main:

        //SALMON_MAP     (trimmedReadPairs, params.salmon_index)

        STAR_GENCODE   (trimmedReadPairs, params.star_index_gencode)

        SALMON_ALIGN   (STAR_GENCODE.out.id, STAR_GENCODE.out.transcriptomeBAM, params.transcriptome_ref)
        
        MULTI_QC       (STAR_GENCODE.out.starLog
                        .collect())

        FILE_CLEANUP   (STAR_GENCODE.out.transcriptomeBAM, STAR_GENCODE.out.starLog, STAR_GENCODE.out.id,MULTI_QC.out.multiqc,SALMON_ALIGN.out.salmonAlignOutput)

}

