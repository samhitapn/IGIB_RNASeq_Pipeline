#!/usr/bin/bash

'RUN SCRIPT => ./nextflowRun.sh \
               /lustre/smriti.a/samhita/src/rnaseqCounts_pipeline/runStats \
                optimisingRun_29032025_mapQuant_test \
                /lustre/smriti.a/samhita/src/rnaseqCounts_pipeline/ \
                false \
                dengue80_all_mappingQuant_test \
                /lustre/smriti.a/samhita/src/rnaseqCounts_pipeline/dengue80.config \
                /lustre/smriti.a/samhita/src/rnaseqCounts_pipeline/alignment_quant_Main.nf '

RUN_LOG_PATH=$1
RUN_CUSTOM_NAME=$2
RESUME=$3
RUNNAME=$4
CONFIG=$5
NF_SCRIPT=$6

cd $1
if [ ! -d "$2" ]; then
    mkdir "$2"
fi

#cd $3

#echo "nextflow run testLearn.nf -c rnaSeq.config -with-trace $1/$2/trace.txt -with-timeline $1/$2/timeline.html -with-report $1/$2/report.html > $1/$2/runLog.log"

if $RESUME; then
    NXF_DEBUG=1 nextflow run $6 \
            -c $5 \
            -with-trace $1/$2/trace.txt \
            -with-timeline $1/$2/timeline.html \
            -with-report $1/$2/report.html \
            -resume \
            --runName $4
else
    NXF_DEBUG=1 nextflow run $6 \
            -c $5 \
            -with-trace $1/$2/trace.txt \
            -with-timeline $1/$2/timeline.html \
            -with-report $1/$2/report.html \
            --runName $4
fi
