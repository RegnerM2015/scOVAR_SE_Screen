ovar_HGSOC_scRNA_processed <- readRDS("./ovar_HGSOC_scRNA_processed.rds")


hgsoc <-subset(ovar_HGSOC_scRNA_processed,cell.type == "0-Epithelial cell")

counts <- hgsoc@assays$RNA@data
counts.bulk <- as.data.frame(rowSums(counts))

counts.bulk$gene <- rownames(counts.bulk)



fimo_outs <- list.files(pattern = "SE60_fimo")


for (i in fimo_outs){
  fimo <- read.table(paste0(i,"/fimo.txt"))
  counts.bulk.new <- counts.bulk[counts.bulk$gene %in% fimo$V2,]
  
  colnames(fimo)[2] <- "gene"
  
  test <- merge(fimo,counts.bulk.new,by = "gene")
  colnames(test)[11] <- "Expr"
  
  test <- dplyr::filter(test,V9 <= 0.1)
  test <- dplyr::arrange(test,desc(Expr))
  write.table(test,paste0(i,"_w_expr.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
}


VlnPlot(ovar_HGSOC_scRNA_processed,features = c("SOX4","ATF4","SOX4","YY1","ATF4","TEAD1"),
        ncol = 6,pt.size = 0,idents = "0-Epithelial cell",same.y.lims = T)+
  ggsave("Enhancer_1_2_3_vln_SE60.pdf",width = 8,height = 4)


ovar_HGSOC_scRNA_processed <- readRDS("./ovar_HGSOC_scRNA_processed.rds")


hgsoc <-subset(ovar_HGSOC_scRNA_processed,cell.type == "11-Epithelial cell")

counts <- hgsoc@assays$RNA@data
counts.bulk <- as.data.frame(rowSums(counts))

counts.bulk$gene <- rownames(counts.bulk)



fimo_outs <- list.files(pattern = "SE14_fimo")


for (i in fimo_outs){
  fimo <- read.table(paste0(i,"/fimo.txt"))
  counts.bulk.new <- counts.bulk[counts.bulk$gene %in% fimo$V2,]
  
  colnames(fimo)[2] <- "gene"
  
  test <- merge(fimo,counts.bulk.new,by = "gene")
  colnames(test)[11] <- "Expr"
  
  test <- dplyr::filter(test,V9 <= 0.1)
  test <- dplyr::arrange(test,desc(Expr))
  write.table(test,paste0(i,"_w_expr.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
}


VlnPlot(ovar_HGSOC_scRNA_processed,features = c("KLF6","YY1","ELF3","KLF6","ELF3","KLF6"),
        ncol = 6,pt.size = 0,idents = "11-Epithelial cell",same.y.lims = T)+
  ggsave("Enhancer_1_2_3_vln_SE14.pdf",width = 8,height = 4)


