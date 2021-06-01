# scOVAR_SE_Screen
1) Cell Ranger was run for the initial processing of the two samples (both scRNA-seq and scATAC-seq). Visit https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome for more details on the CellRanger workflow.

# scRNA-seq processing
1) Run ovar_3BAE2L_RNA-4.R to process the first tumor sample
2) Run ovar_3E5CFL_RNA-4.R to process the second tumor sample
3) Run ovar_HGSOC_RNA.R to merge the two scRNA-seq datasets and reprocess

# scATAC-seq processing
1) Run HGSOC_ATAC.R to process the scATAC-seq data for both samples

# Figure making
1) Run HGSOC_samples_SE_analysis.R to generate the figure visuals
2) Run Motif_Analysis_SE60_SE14.sh to perform the motif analysis
3) Run FIMO_TF_rank.R to process the motif analysis results


NOTE: Similar analyses were conducted at https://github.com/RegnerM2015/scENDO_scOVAR_2020. 


