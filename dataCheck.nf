#!/usr/bin/env nextflow
nextflow.enable.dsl=2

Channel.fromFilePairs("/lustre/smriti.a/samhita/dengue_raw_fastq/*_R{1,2}*.fastq.gz").view()