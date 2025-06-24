## TRIMMOMATIC

   java -jar /lustre/smriti.a/samhita/tools_softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
              -threads ${task.cpus} \
              ${reads.join(' ')} \
              ${sample_id}_R1_paired.fq.gz ${sample_id}_R1_unpaired.fq.gz \
              ${sample_id}_R2_paired.fq.gz ${sample_id}_R2_unpaired.fq.gz \
              ILLUMINACLIP:${params.adapters}:2:20:5 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

## SALMON INDEX GENERATION
module load gcc-8.4

grep "^>" <(gunzip -c /lustre/smriti.a/samhita/reference/GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > /lustre/smriti.a/samhita/reference/salmon_index/salmon_decoys.txt

sed -i.bak -e 's/>//g' /lustre/smriti.a/samhita/reference/salmon_index/salmon_decoys.txt

cat /lustre/smriti.a/samhita/reference/gencode.v47.transcripts.fa.gz /lustre/smriti.a/samhita/reference/GRCh38.primary_assembly.genome.fa.gz > /lustre/smriti.a/samhita/reference/salmon_index/salmon_transcripts_index.fa.gz

salmon index -t /lustre/smriti.a/samhita/reference/salmon_index/salmon_transcripts_index.fa.gz -d /lustre/smriti.a/samhita/reference/salmon_index/salmon_decoys.txt -p 12 -i /lustre/smriti.a/samhita/reference/salmon_index/ --gencode


## SALMON MAPPING MODE

    
    salmon quant -i /lustre/smriti.a/samhita/reference/salmon_index/salmon_decoys.txt -p 12 -i /lustre/smriti.a/samhita/reference/salmon_index/ \
                 -l A \
                 -1 ${trimmedReads[0]} \
                 -2 ${trimmedReads[1]} \
                 --validateMappings \
                 --gcBias \
                 --seqBias \
                 -o ${sample_id}_salmon_quant_mappingMode \
                 -p ${task.cpus}