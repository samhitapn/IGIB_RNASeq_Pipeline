#!/bin/bash -uex
module load gcc-8.4

salmon quant -i salmon_index_keepDuplicates                  -l A                  -1 IGIB11301144292n_Nextseq_Negative_R1_paired.fq.gz                  -2 IGIB11301144292n_Nextseq_Negative_R2_paired.fq.gz                  --validateMappings                  --gcBias                  --seqBias                  -o IGIB11301144292n_Nextseq_Negative_salmon_quant_mappingMode_kd                  -p 10
