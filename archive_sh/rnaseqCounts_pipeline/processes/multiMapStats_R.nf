#!/usr/bin/env nextflow

// Process : Transfer all the output symlinks to the necessary  directories
process MULTI_MAP_STATS {
    
    input:
    path BAM
    val sample_id
    
    publishDir "${params.outDir}/multimappingStats", mode: 'move'
    
    script: 
    
    """
    source /home/smriti.a/anaconda3/etc/profile.d/conda.sh
    conda activate R_m6a_all

    Rscript /lustre/smriti.a/samhita/src/multimapAnalysis.r ${BAM}

    conda deactivate
    """
}

 /*  find . -type l | while read symlink; do
        realfile=$(readlink -f "$symlink")  
        echo $realfile
        symlink_dir=$(dirname "$symlink")
        #echo $symlink_dir
        if [ -d $realfile ]; then
            echo $symlink
            echo $symlink_dir
            rm $symlink
            mv "$realfile" "$symlink_dir/"
        fi
        mv "$realfile" "$symlink_dir"           
    done*/

    /*find ${params.outDir} -type l | while read symlink; do
        realfile=\$(readlink -f "\$symlink")  
        #echo \$realfile
        symlink_dir=\$(dirname "\$symlink")
        if [ -d \$realfile ]; then
            echo \$symlink
            echo \$symlink_dir
            rm \$symlink
            mv "\$realfile" "\$symlink_dir/"
        fi
        mv "\$realfile" "\$symlink_dir/" 
                  
        rm "\$symlink"                       
    done*/