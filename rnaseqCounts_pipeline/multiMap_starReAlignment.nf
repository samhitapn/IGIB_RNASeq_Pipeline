#!/usr/bin/env nextflow
nextflow.enable.dsl=2
//nextflow.preview.output = true

// CHECK FOR MAIN OUTDIR, IF NOT CREATE
def outDir = file(params.outDir)
if (!outDir.exists()) {
  println "Creating output directory: ${params.outDir}"
  outDir.mkdirs()
}

// DEFINE REQUIRED SAMPLES
def reqSamples = ['IGIB113002320529', 'IGIB113002364921', 'IGIB113001484596',
                  'IGIB1130170839','IGIB1130658717','IGIB113001664199',
                  'IGIB1130259871','IGIB113002393366'] as Set

//CREATING A CHANNEL OF PE FILES
trimmedReadPairs = Channel.fromFilePairs(params.trimmedReadsDir, flat: false)
                          /*.toList()
                          .map { list -> list[86,89,105] } 
                          .flatMap { tuple -> tuple }*/
                          .filter { sampleId, reads ->
                               def prefix = sampleId.tokenize('_')[0]
                               //println "SampleID: ${sampleId}, Prefix: ${prefix}"
                               prefix in reqSamples
                         }
                         //.map { sid, reads -> println "SampleID: ${sid}"; [sid, reads] }
                         .view()
                          
                         
//trimmedReadPairs.buffer( size: 50, skip: 0).view()
//trimmedReadPairs.buffer( size: 62, skip: 50).view()

// INCLUDE THE PROCESSES
include { STAR_GENCODE } from './processes/star_gencode_alignment.nf'

// ALL WORKFLOW -> WF_ALL
workflow {
    
    main:

        STAR_GENCODE   (trimmedReadPairs, params.star_index_gencode)

}
