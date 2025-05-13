cat /lustre/smriti.a/samhita/reference/gencode.v47.annotation.gtf | grep "lncRNA" | awk '{ if ($5 - $4 > 200) print }' | wc -l

cat /lustre/smriti.a/samhita/reference/LNCipedia/lncipedia_5_1_hg38.gtf | awk '{ if ($5 - $4 > 200) print }' | wc -l

awk -v RS='>[^\n]+\n' 'length() >= 200 {printf "%s", prt $0} {prt = RT}' /lustre/smriti.a/samhita/reference/gencode.v47.transcripts.fa

seqkit seq -m 200 /lustre/smriti.a/samhita/reference/LNCipedia/lncipedia_5_2.fasta | seqkit stats
seqkit grep -p "\|lncRNA\|" /lustre/smriti.a/samhita/reference/gencode.v47.transcripts.fa | seqkit seq -m 200 | seqkit stats

awk '/^>/ {p=($0 ~ /lncRNA/)} p' /lustre/smriti.a/samhita/reference/gencode.v47.transcripts.fa | seqkit seq -m 200 | awk '/^>/ {header=$0} /lncRNA/ {print header; print}' | awk 'BEGIN {RS=">"; ORS=""} NR>1 {print ">"$0}' | awk -F'|' 'BEGIN {OFS="\t"; print "Transcript_ID", "Gene_ID", "Transcript_Name", "Gene_Name", "Length", "Biotype"}
           {split($1, a, "|"); print substr($1,2), $2, $5, $6, $7, $8}' > /lustre/smriti.a/samhita/reference/lncRNA_ids.txt


           samtools view /lustre/smriti.a/samhita/results_dengue_80/IGIB1130430586n_Nextseq_Negative/star_gencode/IGIB1130430586n_Nextseq_Negative_Aligned.sortedByCoord.out.bam |
  awk 'BEGIN{OFS="\t"} /^@/ || ($0 ~ /NH:i:[2-9]/ || $0 ~ /NH:i:[1-9][0-9]+/)' | less

  awk '{ for(i=12;i<=NF;i++) { if($i ~ /^NH:i:1/) { split($i,a,":"); print a[3] } }}

  samtools view /lustre/smriti.a/samhita/results_dengue_80/IGIB1130430586n_Nextseq_Negative/star_gencode/IGIB1130430586n_Nextseq_Negative_Aligned.sortedByCoord.out.bam chrM | awk '$0 ~ /NH:i:[2-9]/ || $0 ~ /NH:i:[0-9][0-9]+/'