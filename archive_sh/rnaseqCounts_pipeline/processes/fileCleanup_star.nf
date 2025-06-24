#!/usr/bin/env nextflow

// Process : Transfer all the output symlinks to the necessary  directories
process FILE_CLEANUP {
    
    input:
    path transcriptomeBAM
    path starlog
    val sample_id
    //path multiqc
    path salmon

    publishDir "${params.outDir}/${sample_id}/star_gencode", mode: 'copy'
    
    script: 
    
    """
    REAL_BAM=\$(realpath ${transcriptomeBAM})
    REAL_LOG=\$(realpath ${starlog})

    mv \$REAL_BAM ${params.outDir}/${sample_id}/star_gencode/
    mv \$REAL_LOG ${params.outDir}/${sample_id}/star_gencode/
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