#!/bin/bash -uex
module load gcc-8.4

salmon quant -i salmon_index_keepDuplicates                  -l A                  -1 DEN0084_Aviti_Non-Severe_R1_paired.fq.gz                  -2 DEN0084_Aviti_Non-Severe_R2_paired.fq.gz                  --validateMappings                  --gcBias                  --seqBias                  -o DEN0084_Aviti_Non-Severe_salmon_quant_mappingMode_kd                  -p 10
