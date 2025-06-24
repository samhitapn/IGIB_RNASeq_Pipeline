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
bams = Channel.fromPath(params.starBAM)
                          .map { file -> 
                                def sampleId = file.getParent().getParent().getName() // gets the sample folder name
                                return [sampleId, file]
                              }
                          .take(1)
                              
bams
        .map { id, bams -> id }
        .set { sampleID }
        //.view()

bams
        .map { id, bams -> bams  }
        .set { BAM }
        //.view()

include { MULTI_MAP_STATS } from '/lustre/smriti.a/samhita/src/rnaseqCounts_pipeline/processes/multiMapStats_R.nf'


// ALL WORKFLOW -> WF_ALL
workflow {
    
     main:

        MULTI_MAP_STATS   (BAM,sampleID)

}

