########################
# Matt Regner
# Nov 2021
# SE overlap
########################


# Full processed scRNA-seq data was generated in the following directory:
# /datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/scRNA-seq_Processing/HGSOC_RNA

# Fully processed scATAc-seq data was generated in the following directory:
# /datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/scATAC-seq_Processing/HGSOC_ATAC_March_2021-v3



library(ArchR)
library(GenomicRanges)
set.seed(1234)
addArchRGenome(genome="hg38")
addArchRThreads(threads = 16)

###########################
# Subset archr project and call peaks in the malignant fraction 
###########################
proj <- readRDS("./final_archr_proj_archrGS.rds")

levels(factor(proj$predictedGroup_ArchR))

malignant <- c("0-Epithelial cell","2-Epithelial cell","3-Epithelial cell","7-Epithelial cell","11-Epithelial cell","16-Epithelial cell"
)

idxSample <- BiocGenerics::which(proj$predictedGroup_ArchR %in% malignant)
cellsSample <- proj$cellNames[idxSample]
proj <- proj[cellsSample, ]

# Find path to Macs2
pathToMacs2 <- findMacs2()

# Add peakset according to predicted labels
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "predictedGroup_ArchR",force = T)
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "predictedGroup_ArchR",
  pathToMacs2 = pathToMacs2,
  force=T
)
proj <- addPeakMatrix(proj,force=T)

peaks.gr <- getPeakSet(proj)
saveRDS(peaks.gr,"HGSOC_Malignant_Enriched_Peaks_GRanges.rds")

peaks.mat <- getMatrixFromProject(ArchRProj = proj,useMatrix = "PeakMatrix")
saveRDS(peaks.mat,"HGSOC_Malignant_Enriched_Peak_Matrix.rds")

peaks.gr <- data.frame(chrom = peaks.gr@seqnames,
                       ranges = peaks.gr@ranges)
colnames(peaks.gr) <- c("seq","start","end","width","cluster")

write.table(peaks.gr[,1:3],"PeaksOpenInPrimaryTumors.bed",sep = "\t",col.names =F,row.names = F,quote = T)
