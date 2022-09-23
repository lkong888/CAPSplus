#!/bin/bash


# Specify a job name
#$ -N bedbw

# --- Parameters for the Queue Master ---
# Project name and target queue
#$ -P ludwig.prjc
#$ -q short.qc
#$ -pe shmem 4

# Run the job in the current working directory
#$ -cwd 
#$-j y

## filter by depth
d=5 # depth_cutoff
for i in `ls brain*astair_*.filtered.bedGraph` 
do
    grep ^chr[1-9,X,Y] $i|awk -v d=$d '$6>=d' >${i/.bedGraph/}.chr.bedGraph
done

# tab-seq from GSE46710
## convert hg19 to hg38
ln -s ~/cfo155/cfDNA/012020_cfDNA/resource/salk_human_epi/hg19ToHg38.over.chain 
CrossMap.py bed hg19ToHg38.over.chain \
    <(zcat GSE46710_Ad_Front.hmC_sites_FDR_0.01.txt.gz|awk 'BEGIN{OFS="\t"}{if(NR>1)print $1,$2-1,$2,$8,$9,$10,$11}') >GSE46710_Ad_Front.hmC_sites_FDR_0.01.hg38.bed
bedtools intersect -a <(cat ../resource/GSE46710_Ad_Front.hmC_sites_FDR_0.01.hg38.bed|cut -f9-|sort -k1,1 -k2,2n ) \
    -b <(sort -k1,1 -k2,2n brain_healthy_29Jul2022_S3.md.CpG.filtered.chr.bedGraph ) -sorted -wa -wb >brain_healthy_29Jul2022_S3.md.GSE46710_Ad_Front.txt

