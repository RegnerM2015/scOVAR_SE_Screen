########################
# Matt Regner
# March-December 2021
# Super Enhancer Figure
########################


# Full processed scRNA-seq data was generated in the following directory:
# /datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/scRNA-seq_Processing/HGSOC_RNA

# Fully processed scATAc-seq data was generated in the following directory:
# /datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/scATAC-seq_Processing/HGSOC_ATAC_March_2021-v3


library(Seurat)
library(ArchR)
library(ggplot2)
library(stringr)
library(dplyr)
library(forcats)
set.seed(1234)
# Read in color pals:
cols.df <- read.table("Color_Pals.tsv",sep = "\t")

###################################
# Plot marker genes in scRNA-seq
###################################

# Read in scRNA data (processed)
rna <- readRDS("./ovar_HGSOC_scRNA_processed.rds")

# Manually set order of cell type labels
# Define an order of cluster identities
Idents(rna) <- "seurat_clusters"
my_levels <- c("2","3","0","7","11","15","16","19",
               "1","4","9","10","22",
               "14",
               "6","8","13","20","23",
               "5","18","21","12",
               "17")

# Relevel object@ident
rna@active.ident <- factor(x =rna@active.ident, levels = my_levels)

# DimPlot(rna,label = T,label.size = 6)+NoLegend()+ggsave("HGSOC_RNA.pdf",width = 10,height = 10)
# DimPlot(rna,group.by = "Sample")+NoLegend()+ggsave("HGSOC_Sample.pdf",width = 10,height = 10)

p4 <- FeaturePlot(rna,features = "MUC16",cols = c("ivory3","blue"),max.cutoff = 1)+NoLegend()+NoAxes()+theme(title  = element_text(size = 8))
p5 <- FeaturePlot(rna,features = "GJB3",cols = c("ivory3","blue"),max.cutoff = 1)+NoLegend()+NoAxes()+theme(title  = element_text(size = 8))
p6 <- FeaturePlot(rna,features = "IL18",cols = c("ivory3","blue"),max.cutoff = 1)+NoLegend()+NoAxes()+theme(title  = element_text(size = 8))
p7 <- FeaturePlot(rna,features = "LAMB3",cols =c("ivory3","blue"),max.cutoff = 1)+NoLegend()+NoAxes()+theme(title  = element_text(size = 8))
p8 <- FeaturePlot(rna,features = "SERPINB8",cols = c("ivory3","blue"),max.cutoff = 1)+NoLegend()+NoAxes()+theme(title  = element_text(size = 8))
p9 <- FeaturePlot(rna,features = "PMEPA1",cols=c("ivory3","blue"),max.cutoff = 1)+NoLegend()+NoAxes()+theme(title  = element_text(size = 8))
p10 <- FeaturePlot(rna,features = "RELN",cols=c("ivory3","blue"),max.cutoff = 1)+NoLegend()+NoAxes()+theme(title  = element_text(size = 8))
p11 <- FeaturePlot(rna,features = "ROS1",cols=c("ivory3","blue"),max.cutoff = 1)+NoLegend()+NoAxes()+theme(title  = element_text(size = 8))
p12 <- FeaturePlot(rna,features = "SMAD6",cols=c("ivory3","blue"),max.cutoff = 1)+NoLegend()+NoAxes()+theme(title  = element_text(size = 8))


p4.5 <- VlnPlot(rna,features = "MUC16",pt.size = F)+NoLegend()+coord_flip()
p5.5 <- VlnPlot(rna,features = "GJB3",pt.size = F)+NoLegend()+coord_flip()
p6.5 <- VlnPlot(rna,features = "IL18",pt.size = F)+NoLegend()+coord_flip()
p7.5 <- VlnPlot(rna,features = "LAMB3",pt.size = F)+NoLegend()+coord_flip()
p8.5 <- VlnPlot(rna,features = "SERPINB8",pt.size = F)+NoLegend()+coord_flip()
p9.5 <- VlnPlot(rna,features = "PMEPA1",pt.size = F)+NoLegend()+coord_flip()
p10.5 <- VlnPlot(rna,features = "RELN",pt.size = F)+NoLegend()+coord_flip()
p11.5 <- VlnPlot(rna,features = "ROS1",pt.size = F)+NoLegend()+coord_flip()
p12.5 <- VlnPlot(rna,features = "SMAD6",pt.size = F)+NoLegend()+coord_flip()
kruskal.test(p12.5$data$SMAD6 ~ p12.5$data$ident)# Test for signficance 
FeaturePlot(rna,features = "MUC16",cols = c("ivory3","blue"),max.cutoff = 1)+ggsave("Markers-CA125.pdf",width =6,height = 4)
CombinePlots(list(p4,p5,p6,p7,p8,p9,p10,p11,p12),ncol = 2)+ggsave("Markers-feats.pdf",width = 3.5,height = 7.5)




source("stacked_violin.R")
Idents(rna) <- "seurat_clusters"
remove <- c("15","20","19","23","22")#Clusters that did not transfer to ATAC
rna.sub <- subset(x = rna, idents =setdiff(unique(rna$seurat_clusters),remove))
my_levels <- c("2","3","0","7","11","16",
                   "1","4","9","10",
                   "14",
                   "6","8","13",
                   "5","18","21","12",
                   "17")

# Relevel object@ident
rna.sub@active.ident <- factor(x =rna.sub@active.ident, levels = my_levels)
# StackedVlnPlot(rna.sub,features = c("MUC16","GJB3","IL18","LAMB3","SERPINB8","PMEPA1"))+ggsave("Stacked_Violin_SE60.pdf",width = 8,height = 16)
# StackedVlnPlot(rna.sub,features = c("MUC16","SMAD6"))+ggsave("Stacked_Violin_SE14.pdf",width = 8,height =5)
VlnPlot(rna.sub,features = "MUC16")+ggsave("VLN_row_labels.pdf")
StackedVlnPlot(rna.sub,features = c("MUC16","RAE1","EPHA2"))+ggsave("Stacked_Violin_SE60_SE14_combined.pdf",width = 9,height = 8)

# FindMarkers from Seurat
markers <- FindMarkers(rna.sub,ident.1=c(2,3,0,7,11,16),
                       ident.2=c(1, 4, 9, 10, 14, 6, 8, 13, 5, 18, 21, 12, 17),
                       only.pos = T,
                       logfc.threshold = 0)
markers$gene <- rownames(markers)

muc16.test <- VlnPlot(rna.sub,features = "MUC16")
res <- kruskal.test(data=muc16.test$data,MUC16 ~ ident)
print(res)
print(res$p.value)
head(markers[markers$gene == "MUC16",])

rae1.test <- VlnPlot(rna.sub,features = "RAE1")
res <- kruskal.test(data=rae1.test$data,RAE1 ~ ident)
print(res)
print(res$p.value)
head(markers[markers$gene == "RAE1",])

epha2.test <- VlnPlot(rna.sub,features = "EPHA2")
res <- kruskal.test(data=epha2.test$data,EPHA2 ~ ident)
print(res)
print(res$p.value)
head(markers[markers$gene == "EPHA2",])
##########################
# Plot RNA and ATAC UMAPs
##########################
rna <- readRDS("./ovar_HGSOC_scRNA_processed.rds")
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna@meta.data)))
rna.df$Sample <- rna$Sample

rna.sample.plot <-ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = Sample))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = c("dodgerblue3","darkorange2"))+
  guides(colour = guide_legend(override.aes = list(size=6)))
rna.sample.plot +ggsave("Sample_RNA.pdf",width = 8,height = 6)

# ATAC UMAP second:
atac <- readRDS("./proj_LSI_GeneScores_Annotations_Int.rds")
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "Sample",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
atac.df$Sample <- gsub(".*-","",atac.df$color)

ggplot(atac.df,aes_string(x = "x",y="y",color = "Sample"))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values =c("dodgerblue3","darkorange2"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggsave(paste0("Sample_ATAC.pdf"),width = 8,height = 6)
###########################################################




# Plot cell type UMAPs for RNA/ATAC
###########################################################
# RNA UMAP first
rna <- readRDS("./ovar_HGSOC_scRNA_processed.rds")
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna@meta.data)))
rna.df$cell.type <- rna$cell.type

rna.df$cell.type <- factor(rna.df$cell.type,levels = c("2-Epithelial cell","3-Epithelial cell","0-Epithelial cell","7-Epithelial cell","11-Epithelial cell","15-Empty/Epithelial cell","16-Epithelial cell","19-Epithelial cell",
                                                       "1-Fibroblast","4-Fibroblast","9-Fibroblast","10-Fibroblast","22-Fibroblast",
                                                       "14-Endothelial cell",
                                                       "5-T cell","18-T cell","21-T cell","12-NK cell",
                                                       "6-Macrophage","8-Macrophage","13-Macrophage","20-Macrophage","23-Macrophage",
                                                       "17-B cell"
))

cols <- c(cols.df$orange.pal[ seq(6,by=2, len=8)],cols.df$green.pal[seq(4,by=3, len=5)],cols.df$pink.pal[12],cols.df$blue.pal[seq(8,by=3,len=4)],cols.df$grey.pal[seq(6,by=3, len=5)],"black")


rna.sample.plot <-ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color =cell.type))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=3)))
rna.sample.plot +ggsave("Cell_type_RNA.pdf",width = 8,height = 6)

# ATAC UMAP second:
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup_ArchR",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
atac.df$cell.type <- sub(".*?-","",atac.df$color)


atac.df$cell.type <- factor(atac.df$cell.type,levels = c("2-Epithelial cell","3-Epithelial cell","0-Epithelial cell","7-Epithelial cell","11-Epithelial cell","15-Empty_Epithelial cell","16-Epithelial cell","19-Epithelial cell",
                                                         "1-Fibroblast","4-Fibroblast","9-Fibroblast","10-Fibroblast","22-Fibroblast",
                                                         "14-Endothelial cell",
                                                         "5-T cell","18-T cell","21-T cell","12-NK cell",
                                                         "6-Macrophage","8-Macrophage","13-Macrophage","20-Macrophage","23-Macrophage",
                                                         "17-B cell"
))

cols <- c(cols.df$orange.pal[ seq(6,by=2, len=8)],cols.df$green.pal[seq(4,by=3, len=5)],cols.df$pink.pal[12],cols.df$blue.pal[seq(8,by=3,len=4)],cols.df$grey.pal[seq(6,by=3, len=5)],"black")

ggplot(atac.df,aes_string(x = "x",y="y",color = "cell.type"))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ggsave(paste0("Cell_type_ATAC.pdf"),width = 8,height = 6)
###########################################################




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



p <- plotBrowserTrack(
  ArchRProj = atac, 
  groupBy = "predictedGroup_ArchR", 
  region = GRanges("chr20:53730512-53757381"),
  features = GRangesList(TrackA = makeGRangesFromDataFrame(encode),
                         TrackB = makeGRangesFromDataFrame(encode.distal),
                         TrackC = makeGRangesFromDataFrame(encode.epi),
                         TrackD = makeGRangesFromDataFrame(common.snps)),
  loops = NULL,
  upstream = 0,
  downstream = 0
)

dev.off()
pdf("browser_SE60-Extend.pdf",width = 8,height = 14)
grid::grid.draw(p)
dev.off()

p <- plotBrowserTrack(
  ArchRProj = atac, 
  groupBy = "predictedGroup_ArchR", 
  region = GRanges("chr20:53735512-53752381"),
  features = GRangesList(TrackA = makeGRangesFromDataFrame(encode),
                         TrackB = makeGRangesFromDataFrame(encode.distal),
                         TrackC = makeGRangesFromDataFrame(encode.epi),
                         TrackD = makeGRangesFromDataFrame(common.snps)),
  loops = NULL,
  upstream = 0,
  downstream = 0
)

dev.off()
pdf("browser_SE60-Exact.pdf",width = 8,height = 14)
grid::grid.draw(p)
dev.off()


p <- plotBrowserTrack(
  ArchRProj = atac, 
  groupBy = "predictedGroup_ArchR", 
  region = GRanges("chr1:16167176-16190191"),
  loops = NULL,
  features = GRangesList(TrackA = makeGRangesFromDataFrame(encode),
                         TrackB = makeGRangesFromDataFrame(encode.distal),
                         TrackC = makeGRangesFromDataFrame(encode.epi),
                         TrackD = makeGRangesFromDataFrame(common.snps)),
  upstream = 0,
  downstream = 0
)


dev.off()
pdf("browser_SE14-Exact.pdf",width = 8,height = 14)
grid::grid.draw(p)
dev.off()


p <- plotBrowserTrack(
  ArchRProj = atac, 
  groupBy = "predictedGroup_ArchR", 
  region = GRanges("chr1:16162176-16195191
"),
  loops = NULL,
  features = GRangesList(TrackA = makeGRangesFromDataFrame(encode),
                         TrackB = makeGRangesFromDataFrame(encode.distal),
                         TrackC = makeGRangesFromDataFrame(encode.epi),
                         TrackD = makeGRangesFromDataFrame(common.snps)),
  upstream = 0,
  downstream = 0
)


dev.off()
pdf("browser_SE14-Extend.pdf",width = 8,height = 14)
grid::grid.draw(p)
dev.off()

# Find differentially accessible (cancer-enriched peaks)
atac$Cancer.Group <- ifelse(atac$predictedGroup_ArchR %in% c("3_0-Epithelial cell","1_2-Epithelial cell",
                                                             "2_3-Epithelial cell","4_7-Epithelial cell",
                                                             "5_11-Epithelial cell","6_16-Epithelial cell"),"Cancer","Normal")
set.seed(1234)
atac <- addGroupCoverages(ArchRProj = atac, groupBy = "predictedGroup_ArchR",force=T)

atac <- addReproduciblePeakSet(
ArchRProj = atac,
    groupBy = "predictedGroup_ArchR",
    pathToMacs2 = pathToMacs2,force = T
 )
atac <- addPeakMatrix(atac,force = T)

markersPeaks <- getMarkerFeatures(
  ArchRProj = atac,
  useMatrix = "PeakMatrix",
  groupBy = "Cancer.Group",
  useGroups = "Cancer",
  bgdGroups = "Normal",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks,"markersPeaks.rds")

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.25")
cancer <- as.data.frame(markerList$`Cancer`)
cancer$name <- paste0(cancer$seqnames,":",cancer$start,"-",cancer$end)

# SE60 sig peaks
idx <- grep("chr20:537",cancer$name)
cancer.df <- cancer[idx,]
saveRDS(cancer.df,"SigPeaksNearSE60.rds")
cancer.gr <- makeGRangesFromDataFrame(cancer.df)

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
pdf("browser_SE60-Extend-SigPeaks.pdf",width = 8,height = 14)
grid::grid.draw(p)
dev.off()

# SE14 sig peaks
idx <- grep("chr1:161",cancer$name)
cancer.df <- cancer[idx,]
saveRDS(cancer.df,"SigPeaksNearSE14.rds")
cancer.gr <- makeGRangesFromDataFrame(cancer.df)


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
pdf("browser_SE14-Extend-SigPeaks.pdf",width = 8,height = 14)
grid::grid.draw(p)
dev.off()

writeLines(capture.output(sessionInfo()), "sessionInfo-HGSOC_Differential_Genes_And_Peaks.txt")

