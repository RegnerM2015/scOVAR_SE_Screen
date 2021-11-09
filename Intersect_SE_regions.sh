#!/bin/bash
#SBATCH --job-name intersect
#SBATCH --cpus-per-task 8
#SBATCH -c 8
#SBATCH --mem 32g
#SBATCH --partition allnodes

#export R_LIBS_USER="/home/regnerm/R/x86_64-pc-linux-gnu-library/4.0:${R_LIBS_USER}"
source /home/regnerm/anaconda3/etc/profile.d/conda.sh
conda activate scENDO_scOVAR_backup

awk 'BEGIN{FS=OFS="\t"} {gsub(/"/,"\"\"",$1); $1="\"" $1 "\""} 1' hglft_hg38_OVCAR3_BRD4H3K27ac_SE_regions_MikeK.bed > input.bed

bedtools intersect -wa -a input.bed -b PeaksOpenInPrimaryTumors.bed | uniq > output.bed

awk 'BEGIN { FS=","; OFS="\t" } { gsub("\"", "") } { $1=$1 } 1' output.bed > hglft_hg38_OVCAR3_BRD4H3K27ac_SE_regions_MikeK-OpenInPrimaryTumors.bed
