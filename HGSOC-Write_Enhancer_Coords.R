########################
# Matt Regner
# March-December 2021
# Write enhancer coords
########################

library(Seurat)
library(ArchR)
library(ggplot2)
library(stringr)
library(dplyr)
library(forcats)
set.seed(1234)

# Read in atac data:


###########################
# Plot SE browser tracks
###########################
atac <- readRDS("./final_archr_proj_archrGS.rds")


encode <- read.table("./GRCh38-ccREs.bed")
encode <- encode[,1:3]
colnames(encode) <- c("seqname","start","end")
encode.distal <- read.table("./GRCh38-ccREs.dELS.bed")
encode.distal <- encode.distal[,1:3]
colnames(encode.distal)[1:3] <- c("seqname","start","end")
encode.epi <- read.table("./ENCODE_Epithelium_DNase.bed",header = T)
encode.epi <- encode.epi[,1:3]
colnames(encode.epi)[1:3] <- c("seqname","start","end")
common.snps <- read.table("./dbSnp153Common_hg38.bed",sep = "\t")
common.snps <- common.snps[,1:3]
colnames(common.snps)[1:3] <- c("seqname","start","end")
common.snps$end <- common.snps$end+50
clinvar.snps <- read.table("./dbSnp153ClinVar_hg38.bed",sep="\t")
clinvar.snps <- clinvar.snps[,1:3]
colnames(clinvar.snps)[1:3] <- c("seqname","start","end")
clinvar.snps$end <- clinvar.snps$end+50


my_levels <- c("2-Epithelial cell","3-Epithelial cell","0-Epithelial cell","7-Epithelial cell","11-Epithelial cell","16-Epithelial cell",
               "1-Fibroblast","4-Fibroblast","9-Fibroblast","10-Fibroblast",
               "14-Endothelial cell",
               "5-T cell","18-T cell","21-T cell","12-NK cell",
               "6-Macrophage","8-Macrophage","13-Macrophage",
               "17-B cell")

atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "2-Epithelial cell",replacement = "1_2-Epithelial cell")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "3-Epithelial cell",replacement = "2_3-Epithelial cell")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "0-Epithelial cell",replacement = "3_0-Epithelial cell")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "7-Epithelial cell",replacement = "4_7-Epithelial cell")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "11-Epithelial cell",replacement = "5_11-Epithelial cell")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "16-Epithelial cell",replacement = "6_16-Epithelial cell")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "1-Fibroblast",replacement = "7_1-Fibroblast")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "4-Fibroblast",replacement = "8_4-Fibroblast")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "9-Fibroblas",replacement = "9_9-Fibroblas")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "10-Fibroblast",replacement = "10_10-Fibroblast")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "14-Endothelial cell",replacement = "11_14-Endothelial cell")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "6-Macrophage",replacement = "12_6-Macrophage")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "8-Macrophage",replacement = "13_8-Macrophage")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "13-Macrophage",replacement = "14_13-Macrophage")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "5-T cell",replacement = "15_5-T cell")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "18-T cell",replacement = "16_18-T cell")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "21-T cell",replacement = "17_21-T cell")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "12-NK cell",replacement = "18_12-NK cell")
atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = "17-B cell",replacement = "19_17-B cell")

# Read in se60 coords
se60.enh <- readRDS("SigPeaksNearSE60.rds")
se60.enh$name <- as.character(se60.enh$name)
idx.1 <- grep("chr20:5373",se60.enh$name)
idx.2 <- grep("chr20:5374",se60.enh$name)
idx.3 <- grep("chr20:5375",se60.enh$name)
idx <- c(idx.1,idx.2,idx.3)
se60.enh <- se60.enh[idx,]
se60.enh <- dplyr::arrange(se60.enh,start)

idx <- c(1,2,3,4)
se60.enh <- se60.enh[idx,]
se60.enh <- se60.enh[,c(1,3,4)]

# Combine two peaks
se60.enh$end[2] <- se60.enh$end[3]
se60.enh <- se60.enh[-3,]

cancer.gr <- makeGRangesFromDataFrame(se60.enh)

p <- plotBrowserTrack(
  ArchRProj = atac, 
  groupBy = "predictedGroup_ArchR", 
  region = GRanges("chr20:53730512-53757381"),
  features = GRangesList(TrackA = makeGRangesFromDataFrame(encode),
                         TrackB = cancer.gr),
  loops = NULL,
  upstream = 0,
  downstream = 0
)

dev.off()
pdf("browser_SE60-Extend-SigPeaks-Select.pdf",width = 8,height = 14)
grid::grid.draw(p)
dev.off()

# write SE60 enhancer coords to bed files
se60.enh.1 <- se60.enh[1,]
write.table(se60.enh.1,"Enhancer_1_SE60.bed",row.names = F,col.names = F,quote = F,sep = "\t")
se60.enh.2 <- se60.enh[2,]
write.table(se60.enh.2,"Enhancer_2_SE60.bed",row.names = F,col.names = F,quote = F,sep = "\t")
se60.enh.3 <- se60.enh[3,]
write.table(se60.enh.3,"Enhancer_3_SE60.bed",row.names = F,col.names = F,quote = F,sep = "\t")

# Read in se14 coords
se14.enh <- readRDS("SigPeaksNearSE14.rds")
se14.enh$name <- as.character(se14.enh$name)
idx.1 <- grep("chr1:1616",se14.enh$name)
idx.2 <- grep("chr1:1617",se14.enh$name)
idx.3 <- grep("chr1:1618",se14.enh$name)
idx.4 <- grep("chr1:1619",se14.enh$name)
idx <- c(idx.1,idx.2,idx.3,idx.4)
se14.enh <- se14.enh[idx,]
se14.enh <- dplyr::arrange(se14.enh,start)
se14.enh <- se14.enh[c(2:8),]

idx <- c(2,3,4,7)
se14.enh <- se14.enh[idx,]
se14.enh <- se14.enh[,c(1,3,4)]

# Combine first two peaks
se14.enh$end[1] <- se14.enh$end[2]
se14.enh <- se14.enh[-2,]

cancer.gr <- makeGRangesFromDataFrame(se14.enh)

p <- plotBrowserTrack(
  ArchRProj = atac, 
  groupBy = "predictedGroup_ArchR", 
  region = GRanges("chr1:16162176-16195191"),
  loops = NULL,
  features = GRangesList(TrackA = makeGRangesFromDataFrame(encode),
                         TrackB = cancer.gr),
  upstream = 0,
  downstream = 0
)

dev.off()
pdf("browser_SE14-Extend-SigPeaks-Select.pdf",width = 8,height = 14)
grid::grid.draw(p)
dev.off()

# write SE14 enhancer coords to bed files
se14.enh.1 <- se14.enh[1,]
write.table(se14.enh.1,"Enhancer_1_SE14.bed",row.names = F,col.names = F,quote = F,sep = "\t")
se14.enh.2 <- se14.enh[2,]
write.table(se14.enh.2,"Enhancer_2_SE14.bed",row.names = F,col.names = F,quote = F,sep = "\t")
se14.enh.3 <- se14.enh[3,]
write.table(se14.enh.3,"Enhancer_3_SE14.bed",row.names = F,col.names = F,quote = F,sep = "\t")


writeLines(capture.output(sessionInfo()), "sessionInfo-HGSOC-Write_Enhancer_Coords.txt")

