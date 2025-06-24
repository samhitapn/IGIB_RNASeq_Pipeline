library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(reshape2)

bams <- fread("/lustre/smriti.a/samhita/results_dengue_112/multMap_1000_files.txt",
             stringsAsFactors = FALSE,
             header = FALSE)

reads <- readGAlignments(bams[8,V1], param=ScanBamParam(tag=c("AS", "nM", "NH","HI"), what = c("qname","flag")))
multi_reads <- reads[mcols(reads)$NH > 1]

reads_1 <- readGAlignments(bams[1,V1], param=ScanBamParam(tag=c("AS", "nM", "NH","HI"), what = c("qname","flag")))
multi_reads_1 <- reads[mcols(reads)$NH > 1]

multReads_all <- lapply(bams[,V1], function(x){
    sName <- strsplit(x,
                     "/")[[1]][6]
    print(sName)
    print(x)
    getMMReads(x)
    #reads <- readGAlignments(x, param=ScanBamParam(tag=c("AS", "nM", "NH","HI"), what = c("qname","flag")))
    #multi_reads <- reads[mcols(reads)$NH > 1]
    mcols(multi_reads)$Sample <- sName
    print(paste0("Done",sName))
    return(multi_reads)
})

multReads_all_combined <- lapply(multReads_all,function(x){
    data <- as.data.table(x)
    return(data)
})

multReads_all_combined <- rbindlist(multReads_all_combined)

dcast.data.table(multReads_all_combined[,.(count = .N),by = .(Sample,NH)],
                Sample ~ NH)

dcast.data.table(multReads_all_combined[,.(max = max(NH)),by = .(Sample, seqnames)],
                 Sample ~ seqnames)

##### Distribution of the NH per sample

##### Distribution with chr1 to chrM and scaffolds

#####

