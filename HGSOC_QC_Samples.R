########################
# Matt Regner
# Nov 2021
# QC histograms
########################


# Full processed scRNA-seq data was generated in the following directory:
# /datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/scRNA-seq_Processing/HGSOC_RNA

# Fully processed scATAc-seq data was generated in the following directory:
# /datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/scATAC-seq_Processing/HGSOC_ATAC_March_2021-v3


library(Seurat)
library(ArchR)
library(ggplot2)


# Read in scRNA data (processed)
rna <- readRDS("./ovar_HGSOC_scRNA_processed.rds")

rna.df <- rna@meta.data

p1 <- ggplot(rna.df, aes(x=log2(nCount_RNA))) + 
  geom_histogram(binwidth=0.1,color="gray20", fill="gray75")+
  theme_bw()+
  geom_vline(xintercept = log2(500),linetype="dashed",col="darkred")+
  geom_vline(xintercept = log2(median(rna.df$nCount_RNA)),linetype="dashed",col="black")+
  scale_y_continuous(limits=c(0,500))

# Read in scATAC data (processed)
atac <- readRDS("./final_archr_proj_archrGS.rds")

atac.df <- as.data.frame(atac@cellColData)

p2 <- ggplot(atac.df, aes(x=log2(nFrags))) + 
  geom_histogram(binwidth=0.1,color="gray20", fill="gray75")+
  theme_bw()+
  geom_vline(xintercept = log2(1000),linetype="dashed",col="darkred")+
  geom_vline(xintercept = log2(median(atac.df$nFrags)),linetype="dashed",col="black")+
  scale_y_continuous(limits=c(0,500))

p1+p2
ggsave("HGSOC_scRNA_scATAC_QC_Histograms.pdf",width = 7,height = 4)
