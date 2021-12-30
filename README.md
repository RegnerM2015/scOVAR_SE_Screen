# scOVAR_SE_Screen
1) Cell Ranger was run for the initial processing of the two samples (both scRNA-seq and scATAC-seq). Visit https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome for more details on the CellRanger workflow.

# scRNA-seq processing
1) Run ovar_3BAE2L_RNA-4.R to process the first tumor sample
2) Run ovar_3E5CFL_RNA-4.R to process the second tumor sample
3) Run ovar_HGSOC_RNA.R to merge the two scRNA-seq datasets and reprocess

# scATAC-seq processing
1) Run HGSOC_ATAC.R to process the scATAC-seq data for both samples

# Figure making
1) Run HGSOC_Differential_Genes_And_Peaks.R to generate the figure visuals
2) Run HGSOC-Write_Enhancer_Coords.R to write cancer-enriched enhancer bed files for input into motif analysis
3) Run Motif_Analysis_SE60_SE14.sh to perform the motif analysis
4) Run FIMO_TF_rank.R to process the motif analysis results

# Quality control histograms
Run HGSOC_QC_Samples.R to generate histrograms of log2(nCount_RNA) in scRNA-seq and log2(unique fragments) in scATAC-seq

# Overlap peaks called in the malignant fraction with Super Enhancer (SE) regions
1) Run HGSOC_SE_Overlap. R to generate a bed file of peaks called in the malignant fraction (subset to clusters: "0-Epithelial cell","2-Epithelial cell","3-Epithelial cell","7-Epithelial cell","11-Epithelial cell",and "16-Epithelial cell")
2) Run Intersect_SE_regions.sh to find SE regions that overlap with peaks called in the malignant fraction

NOTE: Similar analyses were conducted at https://github.com/RegnerM2015/scENDO_scOVAR_2020. 


