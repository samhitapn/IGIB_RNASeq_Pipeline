java -jar /lustre/smriti.a/samhita/tools_softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
/lustre/smriti.a/samhita/samples/IGIB11301413_R1.fastq.gz /lustre/smriti.a/samhita/samples/IGIB11301413_R2.fastq.gz \
/lustre/smriti.a/samhita/results/trimmomatic/IGIB11301413_R1_paired.fq /lustre/smriti.a/samhita/results/trimmomatic/IGIB11301413_R1_unpaired.fq \
/lustre/smriti.a/samhita/results/trimmomatic/IGIB11301413_R2_paired.fq /lustre/smriti.a/samhita/results/trimmomatic/IGIB11301413_R2_unpaired.fq \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

