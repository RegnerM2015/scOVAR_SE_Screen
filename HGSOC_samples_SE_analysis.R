########################
# Matt Regner
# March 2021
# Super Enhancer Figure
########################


# Full processed scRNA-seq data was generated in the following directory:
# /datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/scRNA-seq_Processing/HGSOC_RNA

# Fully processed scATAc-seq data was generated in the following directory:
# /datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/scATAC-seq_Processing/HGSOC_ATAC_March_2021-v3


library(Seurat)
library(ArchR)
library(ggplot2)


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

FeaturePlot(rna,features = "MUC16",cols = c("ivory3","blue"),max.cutoff = 1)+ggsave("Markers-CA125.pdf",width =6,height = 4)
CombinePlots(list(p4,p5,p6,p7,p8,p9,p10,p11,p12),ncol = 2)+ggsave("Markers-feats.pdf",width = 3.5,height = 7.5)




source("stacked_violin.R")
Idents(rna) <- "seurat_clusters"
my_levels <- c("2","3","0","7","11","15","16","19",
                   "1","4","9","10","22",
                   "14",
                   "6","8","13","20","23",
                   "5","18","21","12",
                   "17")

# Relevel object@ident
rna@active.ident <- factor(x =rna@active.ident, levels = my_levels)
StackedVlnPlot(rna,features = c("MUC16","GJB3","IL18","LAMB3","SERPINB8","PMEPA1"))+ggsave("Stacked_Violin_SE60.pdf",width = 8,height = 16)
StackedVlnPlot(rna,features = c("MUC16","SMAD6"))+ggsave("Stacked_Violin_SE14.pdf",width = 8,height =5)
VlnPlot(rna,features = "MUC16")+ggsave("VLN_row_labels.pdf")
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
pdf("browser_SE60.pdf",width = 8,height = 14)
grid::grid.draw(p)
dev.off()

chr1:16167176-16190191
chr1:16162176-16195191
chr1:16167176-16190191
p <- plotBrowserTrack(
  ArchRProj = atac, 
  groupBy = "predictedGroup_ArchR", 
  region = GRanges("chr1:16167176-16190191
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
pdf("browser_SE14.pdf",width = 8,height = 14)
grid::grid.draw(p)
dev.off()






###########################
# Save peak-to-gene links
###########################
atac <- readRDS("./final_archr_proj_archrGS.rds")

# Read in all P2G peaks 
p2g <- plotPeak2GeneHeatmap(atac,corCutOff = 0,
                                   FDRCutOff = 1,
                                   groupBy = "predictedGroup_ArchR",
                                   returnMatrices = T,nPlot = 5000000)# Make Heatmap object with ALL P2Gs!

mat <- p2g$RNA$matrix
colnames(mat) <- make.unique(p2g$RNA$colData$groupBy)
rownames(mat) <- make.unique(p2g$Peak2GeneLinks$gene)


peaks.genes <- p2g$Peak2GeneLinks

  
peaks.genes <- as.data.frame(peaks.genes)

saveRDS(peaks.genes,"All_P2G_Hits_FDR_of_1.rds")


# SE60 chr20:53735512-53752381
idx <- grep("chr20:537",peaks.genes$peak)
peaks.genes.se60 <- peaks.genes[idx,]
peaks.genes.se60 <- dplyr::filter(peaks.genes.se60,FDR < 0.01)
write.table(peaks.genes.se60,"SE60_linked_genes_BH_FDR_0.01.tsv",sep = "\t",row.names = F,col.names = T,quote = F)

# SE14 chr1:16167176-16190191
idx <- grep("chr1:161",peaks.genes$peak)
peaks.genes.se14 <- peaks.genes[idx,]
peaks.genes.se14 <- dplyr::filter(peaks.genes.se14,FDR < 0.01)
write.table(peaks.genes.se14,"SE14_linked_genes_BH_FDR_0.01.tsv",sep = "\t",row.names = F,col.names = T,quote = F)


enh1 <- data.frame(V1="chr20",V2="53735447",V3="53735947")
write.table(enh1,"enhancer_1_SE60.bed",row.names = F,col.names = F,quote = F,sep = "\t")
enh2 <- data.frame(V1="chr20",V2="53738386",V3="53738886")
write.table(enh2,"enhancer_2_SE60.bed",row.names = F,col.names = F,quote = F,sep = "\t")
enh3 <- data.frame(V1="chr20",V2="53748112",V3="53748612")
write.table(enh3,"enhancer_3_SE60.bed",row.names = F,col.names = F,quote = F,sep = "\t")



enh1 <- data.frame(V1="chr1",V2="16168532",V3="16169032")
write.table(enh1,"enhancer_1_SE14.bed",row.names = F,col.names = F,quote = F,sep = "\t")
enh2 <- data.frame(V1="chr1",V2="16173395",V3="16173895")
write.table(enh2,"enhancer_2_SE14.bed",row.names = F,col.names = F,quote = F,sep = "\t")
enh3 <- data.frame(V1="chr1",V2="16188915",V3="16189415")
write.table(enh3,"enhancer_3_SE14.bed",row.names = F,col.names = F,quote = F,sep = "\t")

#######################
# GSVA analysis 
#######################
source("./stacked_violin.R")
library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(stringr)
library(ComplexHeatmap)


# Read in gene sets and store in gmt format:
files <- list.files(pattern = ".txt")

names <- str_remove(files,".txt")

desired_length <- 14 # or whatever length you want
filler <- vector(mode = "list", length = desired_length)
names(filler) <- names
for (i in names(filler)){
  genes <- read.delim(paste0(i,".txt"),header = T)
  
  genes.vector <- genes$GeneName
  
  filler[[i]] <- genes.vector
  
}

gset <- filler
# Set labels of cell types of interest


# Pseudobulk RNA expression
rna <- readRDS("./ovar_HGSOC_scRNA_processed.rds")

# FeaturePlot(rna,features = "Total_CNVs")+ggsave("CNV.pdf")
# VlnPlot(rna,features = "Total_CNVs")+coord_flip()+ggsave("VLN_CNV.pdf",width = 8,height = 12)


rna.counts <- rna@assays$RNA@counts

rna.pseudobulk <- data.frame(rownames(rna))
for (i in levels(factor(rna$cell.type))){
  cells <- rownames(dplyr::filter(rna@meta.data,cell.type == i))
  
  rna.counts.sub <- rna.counts[,colnames(rna.counts) %in% cells]
  
  rna.counts.bulk <- rowSums(as.matrix(rna.counts.sub))
  
  rna.pseudobulk$i <- rna.counts.bulk
  colnames(rna.pseudobulk)[dim(rna.pseudobulk)[2]] <- i
  
}
rownames(rna.pseudobulk) <- rna.pseudobulk[,1]
rna.pseudobulk <- rna.pseudobulk[,-1]
dim(rna.pseudobulk)


res <- GSVA::gsva(as.matrix(rna.pseudobulk),gset.idx.list = gset,method = "gsva",kcdf = "Poisson")
head(res)

pdf("GSVA_Heatmap_raw_poisson.pdf",width = 6,height = 8)
# Make heatmap annotation
Heatmap(scale(res),cluster_rows = T)
dev.off()



# Rlog transform raw counts
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = rna.pseudobulk,
                              colData = data.frame(meta=colnames(rna.pseudobulk)),design = ~ 1)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

rlog.rna <- rlog(dds,blind = T)
dat <- as.matrix(assay(rlog.rna))


res <- GSVA::gsva(dat,gset.idx.list = gset,method = "gsva",kcdf = "Gaussian")
head(res)

pdf("GSVA_Heatmap_Rlog_gaussian.pdf",width = 6,height =10)
# Make heatmap annotation
Heatmap(scale(res),cluster_rows = T)
dev.off()




# Proportion bar chart scRNA-seq:

meta <- rna@meta.data

df <- meta %>% group_by(RNA_snn_res.0.7) %>% count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(rna$cell.type))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels <- c("2","3","0","7","11","15","16","19",
                                                          "1","4","9","10","22",
                                                          "14",
                                                          "6","8","13","20","23",
                                                          "5","18","21","12",
                                                          "17"))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cell.type))) + 
  geom_bar(position="fill", stat="identity")+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  scale_fill_manual(values = sampleColors)+ggsave("Cell_Type_Prop_RNA.pdf",width = 4,height = 8)






# Patient proportion per subcluster in ATAC:
meta <- as.data.frame(atac@cellColData)
meta$predictedGroup_ArchR <- gsub("-.*", "", meta$predictedGroup_ArchR)

df <- meta %>% group_by(predictedGroup_ArchR) %>% count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(atac$predictedGroup_ArchR))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels =c("2","3","0","7","11","16",
                                                     "1","4","9","10",
                                                     "14",
                                                     "6","8","13",
                                                     "5","18","21","12",
                                                     "17"))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cell.type))) + 
  geom_bar(position="fill", stat="identity")+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  scale_fill_manual(values = sampleColors)+ggsave("Cell_Type_Prop_ATAC.pdf",width = 4,height = 8)
