#!/bin/bash


# Specify a job name
#$ -N anno

# --- Parameters for the Queue Master ---
# Project name and target queue
#$ -P ludwig.prjc
#$ -q short.qc
#$ -pe shmem 4

# Run the job in the current working directory
#$ -cwd -j y

# Log locations which are relative to the current
# working directory of the submission
#$ -o astair_caps_filtering.log
#$ -e error_astair_caps_filtering.log


module load R/4.0.3-foss-2020b
Rscript code/caps_filtering.r 2>&1 | tee > logs/caps_analysis_filtering_astair.log

cat astair_output/caps_all.bedGraph | awk '{{OFS="\t"}}{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}}' | sort -k1,1 -k2,2n > astair_output/caps_all_sorted.bedGraph

module load BEDTools/2.30.0-GCC-10.2.0
bedtools map -a resource/mm9_10kb.bed -b astair_output/caps_all_sorted.bedGraph -c 4,5,6,7,8,9,10,11 -o sum -null "NA" > astair_output/caps_all_rawsignals.bed
##chr, start, end, caps_mereg_mod, caps_merge_unmod, caps_pub_mod, caps_pub_unmod, ace_pub_mod, ace_pub_aC, tab_mod, tab_unmod for astair result
