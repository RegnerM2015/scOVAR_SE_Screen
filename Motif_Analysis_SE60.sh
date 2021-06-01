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

# Get sequence of enhancers and promoter
bedtools getfasta -fi ${dir}/genome.fa -fo enhancer_1.fa -bed enhancer_1.bed
bedtools getfasta -fi ${dir}/genome.fa -fo enhancer_2.fa -bed enhancer_2.bed
bedtools getfasta -fi ${dir}/genome.fa -fo enhancer_3.fa -bed enhancer_3.bed

# Run FIMO motif scanning
fimo --bgfile motif-file --oc enhancer_1_fimo ./JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme enhancer_1.fa
fimo --bgfile motif-file --oc enhancer_2_fimo ./JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme enhancer_2.fa
fimo --bgfile motif-file --oc enhancer_3_fimo ./JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.meme enhancer_3.fa


