#!/bin/bash


# Specify a job name
#$ -N anno_genome

# --- Parameters for the Queue Master ---
# Project name and target queue
#$ -P ludwig.prjc
#$ -q short.qc
#$ -pe shmem 4

# Run the job in the current working directory
#$ -cwd -j y

# Log locations which are relative to the current
# working directory of the submission
###$ -o output_ref.log
###$ -e error_ref.log

## sort HMM file first 
dir=/gpfs2/well/ludwig/users/ebu571/CAPS_haiqi
#hmm_bed=${dir}/resource/mESC_cStates_HMM.bed.gz
#hmm_out=${dir}/resource/mESC_cStates_HMM.sorted.bed.gz
#zcat $hmm_bed | sort -k1,1 -k2,2n > $hmm_out
#gzip $hmm_out

## sort caps methy data
caps_bed=${dir}/astair_output/CAPS_mESC_merged_bam.md.q10_mCtoT_CpG.chr_rmblack.bedGraph.gz
caps_out=${dir}/astair_output/CAPS_mESC_merged_bam.md.q10_mCtoT_CpG.chr_rmblack.sorted.bed.gz
zcat $caps_bed | sort -k1,1 -k2,2n | gzip > $caps_out

module load BEDTools/2.30.0-GCC-10.2.0
## annoate caps methy data
bedb=${dir}/resource/mESC_cStates_HMM.sorted.bed.gz
beda=${dir}/astair_output/CAPS_mESC_merged_bam.md.q10_mCtoT_CpG.chr_rmblack.sorted.bed.gz
out=${dir}/astair_output/CAPS_mESC_merged_bam.md.chr.meth.CpG_annotated.bed.gz
bedtools intersect -sorted -wo -a $beda -b $bedb | gzip > $out
