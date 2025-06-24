library(data.table)
library(reshape2)
library(DESeq2)
library(ggplot2)
library(dplyr)

fcFilePaths <- list.files("/lustre/smriti.a/samhita/results_dengue_80/",
                    pattern = "*_featureCounts_herv.txt$", recursive = TRUE, full.names = TRUE)

featureCountsData <- lapply(fcFilePaths, function(f){
    data <- fread(f,skip = 1, header = TRUE, stringsAsFactors = FALSE)
    return(data[,c(1,7)])
})

featureCountsData_all <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), featureCountsData)
fwrite(featureCountsData_all,
        "/lustre/smriti.a/samhita/results_dengue_80/fcCounts_all.txt",
        quote = FALSE, row.names = FALSE,sep = "\t")

names(featureCountsData_all) <- gsub("_Aligned.sortedByCoord.out.bam","",names(featureCountsData_all))
hervNames <- c(featureCountsData_all[,1])
featureCountsData_all_mat <- as.matrix(featureCountsData_all[,-1])
row.names(featureCountsData_all_mat) <- unlist(hervNames)

# deResults <- getDeseqResults(featureCountsData_all_mat,
#                             salmon = FALSE,
#                             condition = " ~ Severity",
#                             metadata = metadata,
#                             contrasts = c("Severity_Negative_vs_Non-Severe",
#                                           "Severity_Negative_vs_Severe",
#                                           "Severity_Severe_vs_Non-Severe"))

getDeseqResults <- function(data, salmon = TRUE, condition, metadata, contrasts, fcCutOff = NA, padjCutOff = NA, nSamples = NA){
  
  # DESEQ
  if(salmon){
    ddsData <- DESeqDataSetFromTximport(data,
                                   colData = metadata,
                                   design = formula(condition))
  }else{
    ddsData <- DESeqDataSetFromMatrix(data,
                                      colData = metadata,
                                      design = formula(condition))
    print("Got DATA FROM MATRIX")                                    
  }
  
  
  contrastsDEResults <- lapply(contrasts, function(x){
    #### cut-off percentage per group
    print(x)
    #subsets <- c(strsplit(x,"_")[[1]][c(1,2)])
    #if(is.na(nSamples)){
    #  nSamples <- round(min(metadata[Severity %in% subsets,][,.(count = .N),by = Severity][,count])/2)
    #}
    #nSamples <- 10
    print(nSamples)
    print("#########")
    
    keep <- rowSums(counts(ddsData) >= 10) >= nSamples #### Understand why 10
    ddsData <- ddsData[keep, ]
    
    # DESEQ2
    dds <- DESeq(ddsData)
      
    # DEG Results
    contrastDefine <- c(strsplit(x,"_")[[1]][-3])
    results <- results(dds,contrast = contrastDefine)
    results <- data.table("GeneName" = row.names(results),
                                    as.data.table(results))
    results <- results[,contrast := x][,`-log10(padj)` := -log10(padj)]                                    
    return(results)                                
  })
  print("DE ALL CONTRASTS DONE")
  res_contrast <- rbindlist(contrastsDEResults)
  

  print("DONE ALL CONTRASTS")                                    
  
  print("CONCATENATED ALL RESULTS")                                    
  
  # Counts -> To return per contrast group or all?
  # How is the normalisation done?
  counts <- counts(ddsData)
  # vst normalised
  vsd <- vst(ddsData, blind=FALSE)
  #normCounts <- assay(vsd)
  # deseq results
  
  return(list("Counts" = counts,
              "VSD" = vsd,
              "DEGs" = res_contrast,
              "DDS_Data" = ddsData))
  #return(contrastResults_de)
}

##### Differential expression analysis - genes (PC) ####
metadata <- data.table("Sample" = names(featureCountsData_all)[2:79])
metadata <- metadata[,c("Sample","Sequencer","Severity") := tstrsplit(Sample,"_")]
metadata <- metadata[,Severity := gsub("-","",Severity)]

deResults <- getDeseqResults(featureCountsData_all_mat,
                            salmon = FALSE,
                            condition = " ~ Severity",
                            metadata = metadata,
                            contrasts = c("Severity_NonSevere_vs_Negative",
                                          "Severity_Severe_vs_Negative",
                                          "Severity_Severe_vs_NonSevere"),
                            nSamples = round(min(metadata[,.(count = .N),by = Severity][,count])/2))

png("/lustre/smriti.a/samhita/results_dengue_80/pca_ggplot.png")
pca <- plotPCA(deResults$VSD, intgroup = c("Severity","Sequencer"),returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
pcaPlot <- ggplot(pca, aes(x = PC1, y = PC2, color = Severity, shape = Sequencer)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data") +
  theme_test()
ggsave("/lustre/smriti.a/samhita/results_dengue_80/pca_ggplot.png",
    pcaPlot)
dev.off()

save.image("/lustre/smriti.a/samhita/results_dengue_80/deseq_hervAnalysis_12June2025.RData")
#featureCountsData_all_filtered <- featureCountsData_all[,c(1,2,13,24,28:33,3)]


########## Analysisng DE Results
deResults_analysis <- deResults$DEGs
deResults_analysis <- deResults_analysis[,DiffExpression := ifelse((log2FoldChange > 1.5 & padj <= 0.05),"Up",
                                                  ifelse((log2FoldChange < -1.5 & padj <= 0.05), "Down","NS"))]

fwrite(deResults_analysis,
        "/lustre/smriti.a/samhita/results_dengue_80/deResults_HERV.txt",
        quote = FALSE, row.names = FALSE,sep = "\t")

fwrite(deResults$Counts,
        "/lustre/smriti.a/samhita/results_dengue_80/counts_HERV.txt",
        quote = FALSE, row.names = FALSE,sep = "\t")

fwrite(assay(deResults$VSD),
        "/lustre/smriti.a/samhita/results_dengue_80/normCounts_HERV.txt",
        quote = FALSE, row.names = FALSE,sep = "\t")



