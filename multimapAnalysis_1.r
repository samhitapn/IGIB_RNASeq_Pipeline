library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(reshape2)

### Get BAM FILE
# Get arguments passed from command line
args <- commandArgs(trailingOnly = TRUE)

# Check if an argument was given
if (length(args) < 1) {
  stop("Usage: Rscript myscript.R <path_to_file>")
}

# Store the first argument (file path)
getMMReads <- function(bamPath){
    reads <- readGAlignments(bamPath, param=ScanBamParam(tag=c("AS", "nM", "NH","HI"), what = c("qname","flag")))
    multi_reads <- reads[mcols(reads)$NH > 1]

    return(multi_reads)
}

getOverlapSummary <- function(overlapDT, refType){

}

getOverlaps <- function(reads, gtf, sName){

    # Get the overlaps (all)
    hits <- findOverlaps(reads, gtf)
    read_dt <- as.data.table(reads[queryHits(hits)])
    gtf_dt <- as.data.table(gtf[subjectHits(hits)])

    # Rename for clarity
    setnames(read_dt, old = c("start", "end","width"), new = c("read_start", "read_end","read_width"))
    setnames(gtf_dt, old = c("start", "end","width"), new = c("ref_start", "ref_end","ref_width"))
   
    # Add read and gene indices
    read_dt[, read_idx := queryHits(hits)]
    gtf_dt[, gene_idx := subjectHits(hits)]

    # Merge and calculate percent overlap
    overlap_dt <- cbind("Sample" = sName, read_dt, gtf_dt)

    # Compute overlap length
    overlap_dt[, overlap_start := pmax(read_start, ref_start)]
    overlap_dt[, overlap_end := pmin(read_end, ref_end)]
    overlap_dt[, overlap_width := pmax(0, overlap_end - overlap_start + 1)]

    # Percent gene covered by read (relative to gene length) and viceversa
    overlap_dt[, percent_ref_covered := round(100 * overlap_width / ref_width, 2)]
    overlap_dt[, percent_read_covered := round(100 * overlap_width / read_width, 2)]

    # Manual assignment of within categories
    overlap_dt[,overlapType := ifelse(read_start == ref_start & read_end == ref_end,"Equal",
                                      ifelse(read_start >= ref_start & read_end <= ref_end,"Within",
                                      ifelse(read_start <= ref_start & read_end >= ref_end,"Contains",
                                      ifelse(read_start == ref_start,"StartMatch",
                                      ifelse(read_end == ref_end,"EndMatch","Unknown/Others")))))]
    # overlap_dt[read_start == ref_start & read_end == ref_end] <- "equal"
    # overlap_dt[read_start >= ref_start & read_end <= ref_end & overlap_type == ""] <- "within"
    # overlap_dt[read_start <= ref_start & read_end >= ref_end & overlap_type == ""] <- "contains"
    # overlap_dt[read_start == ref_start & overlap_type == ""] <- "start"
    # overlap_dt[read_end == ref_end & overlap_type == ""] <- "end"
    # overlap_dt[overlap_type == ""] <- "partial"  # fallback

    return(overlap_dt)
}


getMultimapped_typeStats <- function(bam, gene, herv, repeats){
    mmReads <- getMMReads(bam)
    sName <- strsplit(bam,
                    "/")[[1]][6]
    overlaps_gene <- getOverlaps(mmReads, gene, sName)
    overlaps_herv <- getOverlaps(mmReads, herv, sName)
    overlaps_repeats <- getOverlaps(mmReads, repeats, sName)
    
}

bam <- args[1]

####### STATS FROM GTF
gtf <- import("/lustre/smriti.a/samhita/reference/gencode.v47.annotation.gtf")
gene <- gtf[mcols(gtf)$type == "gene"]
herv <- import("/lustre/smriti.a/samhita/reference/hervd.gtf")
#repeats <- import("/lustre/smriti.a/samhita/reference/hg38.fa.out.gz")

##### REpeatMasker Data
rmskData <- fread("/lustre/smriti.a/samhita/reference/hg38.fa.out", skip = 3, fill = TRUE)
colnames(rmskData) <- c("SW_score", "perc_div", "perc_del", "perc_ins",
                    "seqname", "start", "end", "left", "strand",
                    "repeat", "class_family", "r_start", "r_end", "r_left", "ID")
rmskData <- rmskData[,strand := ifelse(strand == "C","-",strand)]
repeats <- makeGRangesFromDataFrame(rmskData,keep.extra.columns=TRUE)          

bam <- "/lustre/smriti.a/samhita/results_dengue_112/IGIB113002393366_Severe_M/star_gencode/IGIB113002393366_Severe_M_50_Aligned.sortedByCoord.out.bam"
lapply(bams, function(bam){
    sName <- strsplit(bam,
                    "/")[[1]][6]
    overlaps_gene <- getOverlaps(mmReads, gene, sName)
    overlaps_herv <- getOverlaps(mmReads, herv, sName)
    overlaps_repeats <- getOverlaps(mmReads, repeats, sName)
})
gene <- getMultimapped_typeStats(bam,gene)
herv <- getMultimapped_typeStats(bam)


################################
getMultimapped_typeStats <- function(bam, gtf, herv, rm){
    ####### Filter the multimapped BAMs alone 

    #bam <- "/lustre/smriti.a/samhita/results_dengue_112/IGIB113002393366_Severe_M/star_gencode/IGIB113002393366_Severe_M_Aligned.sortedByCoord.out.bam"
    #bam <- "/lustre/smriti.a/samhita/results_dengue_112/IGIB113002393366_Severe_M/star_gencode/IGIB113002393366_Severe_M_50_Aligned.sortedByCoord.out.bam"
    sName <- strsplit(bam,
                     "/")[[1]][6]
    reads <- readGAlignments(bam, param=ScanBamParam(tag=c("AS", "nM", "NH","HI"), what = c("qname","flag")))
    multi_reads <- reads[mcols(reads)$NH > 1]
   
    #VH00677:1:AAC2MV7M5:1:2406:42204:51622

    hits <- findOverlaps(multi_reads, gtf)
    read_dt <- as.data.table(multi_reads[queryHits(hits)])
    gene_dt <- as.data.table(gtf[subjectHits(hits)])

    # Rename for clarity
    setnames(read_dt, old = c("start", "end","width"), new = c("read_start", "read_end","read_width"))
    setnames(gene_dt, old = c("start", "end","width"), new = c("ref_start", "ref_end","ref_width"))
   
    # Add read and gene indices
    read_dt[, read_idx := queryHits(hits)]
    gene_dt[, gene_idx := subjectHits(hits)]

    # Step 3: Merge and calculate percent overlap
    overlap_dt <- cbind(read_dt, gene_dt)

    # Compute overlap length
    overlap_dt[, overlap_start := pmax(read_start, ref_start)]
    overlap_dt[, overlap_end := pmin(read_end, ref_end)]
    overlap_dt[, overlap_width := pmax(0, overlap_end - overlap_start + 1)]

    # Percent gene covered by read (relative to gene length)
    overlap_dt[, percent_ref_covered := round(100 * overlap_width / ref_width, 2)]
    overlap_dt[, percent_read_covered := round(100 * overlap_width / read_width, 2)]

    # Final cleanup
    overlap_summary <- overlap_dt[, .(
        sample_name = sName,
        read_name = qname,
        chromosome = seqnames,
        read_start,
        read_end,
        readWidth = qwidth,
        alignWidth = width,
        flag,
        strand,
        NH,
        HI,
        AS,
        nM,
        CIGAR = cigar,
        gene_id,
        category = "Genes",
        type,
        Biotype = gene_type,
        ref_start,
        ref_end,
        overlap_width,
        ref_width,
        read_width,
        percent_ref_covered,
        percent_read_covered
    )]

    ######### HERV
    hits_herv <- findOverlaps(multi_reads, herv)
    read_dt_herv <- as.data.table(multi_reads[queryHits(hits_herv)])
    gene_dt_herv <- as.data.table(herv[subjectHits(hits_herv)])

    # Rename for clarity
    setnames(read_dt_herv, old = c("start", "end","width"), new = c("read_start", "read_end","read_width"))
    setnames(gene_dt_herv, old = c("start", "end","width"), new = c("ref_start", "ref_end","ref_width"))
   

    # Add read and gene indices
    read_dt_herv[, read_idx := queryHits(hits_herv)]
    gene_dt_herv[, gene_idx := subjectHits(hits_herv)]

    # Step 3: Merge and calculate percent overlap
    overlap_dt_herv <- cbind(read_dt_herv, gene_dt_herv)

    # Compute overlap length
    overlap_dt_herv[, overlap_start := pmax(read_start, ref_start)]
    overlap_dt_herv[, overlap_end := pmin(read_end, ref_end)]
    overlap_dt_herv[, overlap_width := pmax(0, overlap_end - overlap_start + 1)]

    # Percent gene covered by read (relative to gene length)
    overlap_dt_herv[, percent_ref_covered := round(100 * overlap_width / ref_width, 2)]
    overlap_dt_herv[, percent_read_covered := round(100 * overlap_width / read_width, 2)]

    # Final cleanup
    overlap_summary_herv <- overlap_dt_herv[, .(
        sample_name = sName,
        read_name = qname,
        read_start,
        read_end,
        readWidth = qwidth,
        alignWidth = width,
        flag,
        strand,
        NH,
        HI,
        AS,
        nM,
        CIGAR = cigar,
        gene_id,
        category = "HERV",
        type,
        Biotype = "HERVd",
        ref_start,
        ref_end,
        overlap_width,
        ref_width,
        read_width,
        percent_ref_covered,
        percent_read_covered
    )]

    ######### REPEAT MASKER
    hits_rm <- findOverlaps(multi_reads, rm)
    read_dt_rm <- as.data.table(multi_reads[queryHits(hits_rm)])
    gene_dt_rm <- as.data.table(rmskData_gr[subjectHits(hits_rm)])

    # Rename for clarity
    setnames(read_dt_rm, old = c("start", "end","width"), new = c("read_start", "read_end","read_width"))
    setnames(gene_dt_rm, old = c("start", "end","width"), new = c("ref_start", "ref_end","ref_width"))
   

    # Add read and gene indices
    read_dt_rm[, read_idx := queryHits(hits_rm)]
    gene_dt_rm[, gene_idx := subjectHits(hits_rm)]

    # Step 3: Merge and calculate percent overlap
    overlap_dt_rm <- cbind(read_dt_rm, gene_dt_rm)

    # Compute overlap length
    overlap_dt_rm[, overlap_start := pmax(read_start, ref_start)]
    overlap_dt_rm[, overlap_end := pmin(read_end, ref_end)]
    overlap_dt_rm[, overlap_width := pmax(0, overlap_end - overlap_start + 1)]

    # Percent gene covered by read (relative to gene length)
    overlap_dt_rm[, percent_ref_covered := round(100 * overlap_width / ref_width, 2)]
    overlap_dt_rm[, percent_read_covered := round(100 * overlap_width / read_width, 2)]

    # Final cleanup
    overlap_summary_rm <- overlap_dt_rm[, .(
        sample_name = sName,
        read_name = qname,
        read_start,
        read_end,
        readWidth = qwidth,
        alignWidth = width,
        flag,
        strand,
        NH,
        HI,
        AS,
        nM,
        CIGAR = cigar,
        gene_id = repeat.,
        category = "Repeats",
        type = "Repeats",
        Biotype = class_family,
        ref_start,
        ref_end,
        overlap_width,
        ref_width,
        read_width,
        percent_ref_covered,
        percent_read_covered
    )]

    fwrite(overlap_summary_herv,
          paste0("/lustre/smriti.a/samhita/results_dengue_112/multimappingStats/",sName,"_herv.tsv"),
          quote = FALSE,
          row.names = FALSE,sep = "\t")

    fwrite(overlap_summary,
          paste0("/lustre/smriti.a/samhita/results_dengue_112/multimappingStats/",sName,"_genes.tsv"),
          quote = FALSE,
          row.names = FALSE,sep = "\t")

    fwrite(overlap_summary_rm,
          paste0("/lustre/smriti.a/samhita/results_dengue_112/multimappingStats/",sName,"_rm.tsv"),
          quote = FALSE,
          row.names = FALSE,sep = "\t")
    ####### return
}


getMultimapped_typeStats(bam, gtf, herv, rmskData_gr)   

############# ANALYSIS ########
testSample <- fread("/lustre/smriti.a/samhita/results_dengue_112/multimappingStats/IGIB113002393366_Severe_M_genes.tsv",
                    header = TRUE,
                    stringsAsFactors = FALSE)
testSample_herv <- fread("/lustre/smriti.a/samhita/results_dengue_112/multimappingStats/IGIB113002393366_Severe_M_herv.tsv",
                    header = TRUE,
                    stringsAsFactors = FALSE)

testSample_rm <- fread("/lustre/smriti.a/samhita/results_dengue_112/multimappingStats/IGIB113002393366_Severe_M_rm.tsv",
                    header = TRUE,
                    stringsAsFactors = FALSE)

# repeatRef <- import("/lustre/smriti.a/samhita/reference/hg38.fa.out.gz")
# bam <- "/lustre/smriti.a/samhita/results_dengue_80/DEN0010_Aviti_Non-Severe/star_gencode/DEN0010_Aviti_Non-Severe_Aligned.sortedByCoord.out.bam"
# mBam <- "/lustre/smriti.a/samhita/results_dengue_80/DEN0010_Aviti_Non-Severe/star_gencode/DEN0010_Aviti_Non-Severe_Aligned.sortedByCoord_multimapped.bam"

# getMultimapped_typeStats(bam,mBam)


# samtools view -h "/lustre/smriti.a/samhita/results_dengue_80/DEN0010_Aviti_Non-Severe/star_gencode/DEN0010_Aviti_Non-Severe_Aligned.sortedByCoord.out.bam" \
#     | awk 'BEGIN{OFS="\t"} /^@/ || ($0 ~ /NH:i:[2-9]/ || $0 ~ /NH:i:[1-9][0-9]+/)' \
#     | samtools view -b -o "/lustre/smriti.a/samhita/results_dengue_80/DEN0010_Aviti_Non-Severe/star_gencode/DEN0010_Aviti_Non-Severe_Aligned.sortedByCoord_multimapped.bam"

# bedtools bamtobed -bed12 -i "/lustre/smriti.a/samhita/results_dengue_80/DEN0010_Aviti_Non-Severe/star_gencode/DEN0010_Aviti_Non-Severe_Aligned.sortedByCoord_multimapped.bam" > "/lustre/smriti.a/samhita/results_dengue_80/DEN0010_Aviti_Non-Severe/star_gencode/DEN0010_Aviti_Non-Severe_multimapped_bed12.bed"

samtools view -h "/lustre/smriti.a/samhita/results_dengue_112/IGIB113002393366_Severe_M/star_gencode/IGIB113002393366_Severe_M_Aligned.sortedByCoord.out.bam" \
    | awk 'BEGIN{OFS="\t"} /^@/ || ($0 ~ /NH:i:[2-9]/ || $0 ~ /NH:i:[1-9][0-9]+/)' \
    | samtools view -b - \
    | bedtools bamtobed -bed12 -i - > "/lustre/smriti.a/samhita/results_dengue_112/IGIB113002393366_Severe_M/star_gencode/IGIB113002393366_Severe_M_Aligned_multimapped.bed"

# samtools view /lustre/smriti.a/samhita/results_dengue_80/DEN0010_Aviti_Non-Severe/star_gencode/DEN0010_Aviti_Non-Severe_Aligned.sortedByCoord.out.bam | awk '{
#   for (i = 12; i <= NF; i++) {
#     if ($i ~ /^NH:i:/) {
#       split($i, a, ":")
#       if (a[3] > 1 && a[3] > max) max = a[3]
#     }
#   }
# } END { print max }'

bedtools intersect -a "/lustre/smriti.a/samhita/results_dengue_112/IGIB113002393366_Severe_M/star_gencode/IGIB113002393366_Severe_M_Aligned_multimapped.bed" -b "/lustre/smriti.a/samhita/reference/gencode.v47.annotation.gtf.bed" -wa -wb > "/lustre/smriti.a/samhita/results_dengue_112/IGIB113002393366_Severe_M/star_gencode/IGIB113002393366_Severe_M_Aligned_multimapped_intersect.bed"