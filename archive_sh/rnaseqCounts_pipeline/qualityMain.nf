#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// CHECK FOR MAIN OUTDIR, IF NOT CREATE
def outDir = file(params.outDir)
if (!outDir.exists()) {
    println "Creating output directory: ${params.outDir}"
    outDir.mkdirs()
}

//CREATING A CHANNEL OF PE FILES
//rawReadPairs = Channel.fromFilePairs(params.readsDir)
dengue_classification = Channel.fromPath(params.dengue)
       .splitCsv( header: ['col0','col1','col2'], skip: 1 )
       .map(row -> tuple(row.col0,"${row.col0}_${row.col1}_${row.col2}"))
       

rawReadPairs = Channel.fromFilePairs(params.readsDir)
        .map{id, reads -> 
            def id_1 = id.split("_")[0]
            return(tuple(id_1,reads))
        }
        .join(dengue_classification)
        .map{ id_1, reads, id_new -> tuple(id_new, reads)}

// INCLUDE THE PROCESSES
include { FASTQC } from './processes/fastqc.nf'
include { FASTQC as FASTQC_POST } from './processes/fastqc.nf'
include { TRIM } from './processes/trimmomatic.nf'
/**include { SALMON_MAP } from './processes/salmon_mappingMode.nf'
include { STAR_GENCODE } from './processes/star_gencode_alignment.nf'
include { SALMON_ALIGN } from './processes/salmon_alignmentMode.nf'*/
include { MULTI_QC } from './processes/multiqc.nf'
include { FILE_TRANSFER } from './processes/fileTransfer.nf'

// ALL WORKFLOW -> WF_ALL
workflow {
    
    FASTQC         (rawReadPairs.map { id, reads -> tuple(id, reads, "preTrim") })  
    
    TRIM           (rawReadPairs)
  
    FASTQC_POST    (TRIM.out.trimOutput.map {id,R1,R2 -> tuple(id, [R1, R2], "postTrim")})

    //SALMON_MAP     (TRIM.out.trimOutput.map {id,R1,R2 -> tuple(id, [R1, R2])}, params.salmon_index)

    //STAR_GENCODE   (TRIM.out.trimOutput.map {id,R1,R2 -> tuple(id, [R1, R2])}, params.star_index_gencode)

    //SALMON_ALIGN   (STAR_GENCODE.out.id, STAR_GENCODE.out.transcriptomeBAM, params.transcriptome_ref)

  
    MULTI_QC       (FASTQC_POST.out.zip
                     .mix(FASTQC_POST.out.html)
                     .collect())

    allOutput = FASTQC.out.html.mix(FASTQC.out.zip,
                     TRIM.out.R1_up,TRIM.out.R2_up,TRIM.out.trimOutput.map{id, R1, R2 -> [R1,R2]},
                     FASTQC_POST.out.html,FASTQC_POST.out.zip,
                     MULTI_QC.out.multiqc)
                     .collect()      
    FILE_TRANSFER   (allOutput)
    
}

