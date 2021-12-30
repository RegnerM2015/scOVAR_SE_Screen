########################
# Matt Regner
# March-December 2021
# Process Motif Results
########################

library(Seurat)
library(ArchR)
library(ggplot2)
library(stringr)
library(dplyr)

ovar_HGSOC_scRNA_processed <- readRDS("./ovar_HGSOC_scRNA_processed.rds")

malignant <- c("0-Epithelial cell","2-Epithelial cell",
"3-Epithelial cell","7-Epithelial cell",
"11-Epithelial cell","16-Epithelial cell")

hgsoc <-subset(ovar_HGSOC_scRNA_processed,cell.type %in% malignant)

counts <- hgsoc@assays$RNA@data
counts.bulk <- as.data.frame(rowSums(counts))

counts.bulk$gene <- rownames(counts.bulk)



fimo_outs <- list.dirs()[grep("SE60_fimo",list.dirs())]


for (i in fimo_outs){
  fimo <- read.table(paste0(i,"/fimo.txt"))
  counts.bulk.new <- counts.bulk[counts.bulk$gene %in% fimo$V2,]
  
  colnames(fimo)[2] <- "gene"
  
  test <- merge(fimo,counts.bulk.new,by = "gene")
  colnames(test)[11] <- "Expr"
  
  test <- dplyr::filter(test,V9 <= 0.1)
  test <- dplyr::arrange(test,desc(Expr))
  write.table(test,paste0(i,"_w_expr-Update.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
}

hgsoc$Cancer <- ifelse(hgsoc$cell.type %in% malignant,"Cancer","Non-Cancer")
Idents(hgsoc) <- "Cancer"
VlnPlot(hgsoc,features = c("ATF4","SOX4","ELF3","YY1","ATF4","JUN"),
        ncol = 6,pt.size = 0,idents = "Cancer",same.y.lims = T)+
  ggsave("Enhancer_1_2_3_vln_SE60-Update.pdf",width = 8,height = 4)


ovar_HGSOC_scRNA_processed <- readRDS("./ovar_HGSOC_scRNA_processed.rds")

malignant <- c("0-Epithelial cell","2-Epithelial cell",
               "3-Epithelial cell","7-Epithelial cell",
               "11-Epithelial cell","16-Epithelial cell")

hgsoc <-subset(ovar_HGSOC_scRNA_processed,cell.type %in% malignant)

counts <- hgsoc@assays$RNA@data
counts.bulk <- as.data.frame(rowSums(counts))

counts.bulk$gene <- rownames(counts.bulk)



fimo_outs <- list.dirs()[grep("SE14_fimo",list.dirs())]


for (i in fimo_outs){
  fimo <- read.table(paste0(i,"/fimo.txt"))
  counts.bulk.new <- counts.bulk[counts.bulk$gene %in% fimo$V2,]
  
  colnames(fimo)[2] <- "gene"
  
  test <- merge(fimo,counts.bulk.new,by = "gene")
  colnames(test)[11] <- "Expr"
  
  test <- dplyr::filter(test,V9 <= 0.1)
  test <- dplyr::arrange(test,desc(Expr))
  write.table(test,paste0(i,"_w_expr-Update.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
}

hgsoc$Cancer <- ifelse(hgsoc$cell.type %in% malignant,"Cancer","Non-Cancer")
Idents(hgsoc) <- "Cancer"
VlnPlot(hgsoc,features = c("ELF3","KLF5","JUND","JUNB","KLF6","ELF3"),
        ncol = 6,pt.size = 0,idents = "Cancer",same.y.lims = T)+
  ggsave("Enhancer_1_2_3_vln_SE14-Update.pdf",width = 8,height = 4)


writeLines(capture.output(sessionInfo()), "sessionInfo-FIMO_TF_Rank.txt")
