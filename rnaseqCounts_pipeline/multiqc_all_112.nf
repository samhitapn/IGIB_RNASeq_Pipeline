#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// CHECK FOR MAIN OUTDIR, IF NOT CREATE
def outDir = file(params.outDir)
if (!outDir.exists()) {
    println "Creating output directory: ${params.outDir}"
    outDir.mkdirs()
}

//CREATING A CHANNEL OF PE FASTQ FILES
fastqcResults_preTrim = Channel.fromPath("/lustre/smriti.a/samhita/results_dengue_112/**/fastqc/*_preTrim/*")
                               .collect()

fastqcResults_postTrim = Channel.fromPath("/lustre/smriti.a/samhita/results_dengue_112/**/fastqc/*_postTrim/*")
                                .collect()

salmonMapResults = Channel.fromPath("/lustre/smriti.a/samhita/results_dengue_112/**/salmon_map/*_salmon_quant_mappingMode/",type: 'dir')
                          .collect()

salmonAlignResults = Channel.fromPath("/lustre/smriti.a/samhita/results_dengue_112/**/salmon_align/*_salmon_quant_alignmentMode/",type: 'dir')
                          .collect()

starResults = Channel.fromPath("/lustre/smriti.a/samhita/results_dengue_112/**/star_gencode/*_1000_Log.final.out")
                     .collect()

//                    fastqcResults_preTrim.concat(salmonMapResults).concat(fastqcResults_postTrim).concat(starResults).collect().view()


// INCLUDE THE PROCESSES
include { MULTI_QC } from '/lustre/smriti.a/samhita/src/rnaseqCounts_pipeline/processes/multiqc_process.nf'


// ALL WORKFLOW -> WF_ALL
workflow {
    
    MULTI_QC       (starResults)
}

