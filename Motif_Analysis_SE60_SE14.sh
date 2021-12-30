#!/bin/bash
#SBATCH --job-name Motif 
#SBATCH --cpus-per-task 8
#SBATCH -c 8
#SBATCH --mem 32g
#SBATCH --partition allnodes

export PATH=/home/regnerm/Software/cellranger-dna-1.1.0:$PATH
export PATH=/home/regnerm/Software/tabix:$PATH
export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-4.12.0:$PATH
export PATH=/share/apps/install-compute/bin/bcftools:$PATH

dir=/datastore/nextgenout5/share/labs/francolab/Data/refdata-cellranger-atac-GRCh38-1.2.0/fasta

# Get sequence of enhancers for SE60
bedtools getfasta -fi ${dir}/genome.fa -fo enhancer_1_SE60.fa -bed Enhancer_1_SE60.bed
bedtools getfasta -fi ${dir}/genome.fa -fo enhancer_2_SE60.fa -bed Enhancer_2_SE60.bed
bedtools getfasta -fi ${dir}/genome.fa -fo enhancer_3_SE60.fa -bed Enhancer_3_SE60.bed

# Run FIMO motif scanning
fimo --bgfile motif-file --oc enhancer_1_SE60_fimo ./JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme enhancer_1_SE60.fa
fimo --bgfile motif-file --oc enhancer_2_SE60_fimo ./JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme enhancer_2_SE60.fa
fimo --bgfile motif-file --oc enhancer_3_SE60_fimo ./JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme enhancer_3_SE60.fa


# Get sequence of enhancers for SE14
bedtools getfasta -fi ${dir}/genome.fa -fo enhancer_1_SE14.fa -bed Enhancer_1_SE14.bed
bedtools getfasta -fi ${dir}/genome.fa -fo enhancer_2_SE14.fa -bed Enhancer_2_SE14.bed
bedtools getfasta -fi ${dir}/genome.fa -fo enhancer_3_SE14.fa -bed Enhancer_3_SE14.bed

# Run FIMO motif scanning
fimo --bgfile motif-file --oc enhancer_1_SE14_fimo ./JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme enhancer_1_SE14.fa
fimo --bgfile motif-file --oc enhancer_2_SE14_fimo ./JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme enhancer_2_SE14.fa
fimo --bgfile motif-file --oc enhancer_3_SE14_fimo ./JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme enhancer_3_SE14.fa

