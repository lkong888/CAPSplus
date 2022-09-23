#!/bin/bash


# Specify a job name
#$ -N phred

# --- Parameters for the Queue Master ---
# Project name and target queue
#$ -P ludwig.prjc
#$ -q short.qc
#$ -pe shmem 4

# Run the job in the current working directory
#$ -cwd -j y

# Log locations which are relative to the current
# working directory of the submission
#$ -o phred_score.log
#$ -e phred_score.log

module load Python/3.7.4-GCCcore-8.3.0

mESC1_r1=/users/ludwig/ebu571/ebu571/CAPS_haiqi/fastq/mESC_1_29Jul2022_S1_R1.fastq.gz
mESC1_r2=/users/ludwig/ebu571/ebu571/CAPS_haiqi/fastq/mESC_1_29Jul2022_S1_R2.fastq.gz
mESC2_r1=/users/ludwig/ebu571/ebu571/CAPS_haiqi/fastq/mESC_2_29Jul2022_S2_R1.fastq.gz
mESC2_r2=/users/ludwig/ebu571/ebu571/CAPS_haiqi/fastq/mESC_2_29Jul2022_S2_R2.fastq.gz

zcat $mESC1_r1 $mESC2_r1 | gzip > /users/ludwig/ebu571/ebu571/CAPS_haiqi/fastq/mESC_R1.fastq.gz
zcat $mESC1_r2 $mESC2_r2 | gzip > /users/ludwig/ebu571/ebu571/CAPS_haiqi/fastq/mESC_R2.fastq.gz
 
/gpfs2/well/ludwig/users/ebu571/Python/asTair_p3.7.4-skylake/bin/astair phred -1 fastq/mESC_R1.fastq.gz -2 fastq/mESC_R2.fastq.gz --plot -d /users/ludwig/ebu571/ebu571/CAPS_haiqi/temp
