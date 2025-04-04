#!/usr/bin/env nextflow

// Define the input file
dengue80_classification = Channel.fromPath("/lustre/smriti.a/samhita/dengue_raw_fastq/dengue_80_Classification.csv")
       .splitCsv( header: ['col0','col1','col2','col3','col4','col5','col6'], skip: 1 )
       .map(row -> tuple(row.col1,"${row.col1}_${row.col6}_${row.col4}"))

Channel.fromFilePairs("/lustre/smriti.a/samhita/dengue_raw_fastq/DEN00*_R{1,2}*.fastq.gz")
       .map{id, reads -> 
            def id_1 = id.split("_")[0]
            return(tuple(id_1,reads))
        }
        .join(dengue80_classification)
        .map{ id_1, reads, id_new -> tuple(id_new, reads)}
        .view()