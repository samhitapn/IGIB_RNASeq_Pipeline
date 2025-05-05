awk '$3 == "exon"' /lustre/smriti.a/samhita/reference/gencode.v47.annotation.gtf | \
  awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, ".", ".", $7}' | sort -k1,1 -k2,2n -k3,3n -k6,6 | uniq > /lustre/smriti.a/samhita/reference/manual_exons.bed


awk '$3 == "gene"' /lustre/smriti.a/samhita/reference/gencode.v47.annotation.gtf | \
  awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, $10, ".", $7}' | sort -k1,1 -k2,2n -k3,3n -k6,6 | uniq > /lustre/smriti.a/samhita/reference/manual_genes.bed


bedtools subtract -a /lustre/smriti.a/samhita/reference/manual_genes.bed -b /lustre/smriti.a/samhita/reference/manual_exons.bed > /lustre/smriti.a/samhita/reference/manual_introns.bed

samtools index /lustre/smriti.a/samhita/results_dengue_80/IGIB1130746713n_Nextseq_Negative/star_gencode/IGIB1130746713n_Nextseq_Negative_Aligned.sortedByCoord.out.bam

bedtools bamtobed -i /lustre/smriti.a/samhita/results_dengue_80/IGIB1130746713n_Nextseq_Negative/star_gencode/IGIB1130746713n_Nextseq_Negative_Aligned.sortedByCoord.out.bam > /lustre/smriti.a/samhita/results_dengue_80/IGIB1130746713n_Nextseq_Negative/star_gencode/IGIB1130746713n_Nextseq_Negative_reads.bed
sort /lustre/smriti.a/samhita/results_dengue_80/IGIB1130746713n_Nextseq_Negative/star_gencode/IGIB1130746713n_Nextseq_Negative_reads.bed | uniq -d
# stranded exon
bedtools intersect -b /lustre/smriti.a/samhita/results_dengue_80/IGIB1130746713n_Nextseq_Negative/star_gencode/IGIB1130746713n_Nextseq_Negative_reads.bed \
-a /lustre/smriti.a/samhita/reference/manual_exons.bed -s -u | wc -l

# unstranded exon
bedtools intersect -a /lustre/smriti.a/samhita/results_dengue_80/IGIB1130746713n_Nextseq_Negative/star_gencode/IGIB1130746713n_Nextseq_Negative_reads.bed \
-b /lustre/smriti.a/samhita/reference/manual_exons.bed -u | wc -l

# stranded intron
bedtools intersect -b /lustre/smriti.a/samhita/results_dengue_80/IGIB1130746713n_Nextseq_Negative/star_gencode/IGIB1130746713n_Nextseq_Negative_reads.bed \
-a /lustre/smriti.a/samhita/reference/manual_introns.bed -s -u | wc -l

# unstranded intron
bedtools intersect -a /lustre/smriti.a/samhita/results_dengue_80/IGIB1130746713n_Nextseq_Negative/star_gencode/IGIB1130746713n_Nextseq_Negative_reads.bed \
-b /lustre/smriti.a/samhita/reference/manual_introns.bed -u | wc -l
