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
alignedBams = Channel.fromPath(params.starBAM)
                          .map { file -> 
                                def sampleId = file.getParent().getParent().getName() // gets the sample folder name
                                return [sampleId, file]
                              }
                          .view()
/*alignedBams
        .map { id, alignedBams -> id }
        .set { sampleID }

alignedBams
        .map { id, alignedBams -> alignedBams  }
        .set { BAM }

include { FEATURE_COUNTS } from '/lustre/smriti.a/samhita/src/rnaseqCounts_pipeline/processes/featureCounts.nf'


// ALL WORKFLOW -> WF_ALL
workflow {
    
     main:

        FEATURE_COUNTS   (sampleID, BAM, params.herv_gtf)

}*/

