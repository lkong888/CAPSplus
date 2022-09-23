#!/bin/bash


# Specify a job name
#$ -N SI

# --- Parameters for the Queue Master ---
# Project name and target queue
#$ -P ludwig.prjc
#$ -q short.qc
#$ -pe shmem 4

# Run the job in the current working directory
#$ -cwd 
#$-j y

# Log locations which are relative to the current
# working directory of the submission
#$ -o downsample_ref.log
#$ -e downsample_ref.log
module load samtools/1.8-gcc5.4.0
caps_plus="/users/ludwig/ebu571/ebu571/CAPS_haiqi/align/CAPS_mESC_merged_bam.md.q10.bam"
chrbam="/users/ludwig/ebu571/ebu571/CAPS_haiqi/align/CAPS_mESC_merged_bam.md.q10.chr.bam"
down="/users/ludwig/ebu571/ebu571/CAPS_haiqi/align/CAPS_mESC_merged_bam.md.q10.chr_0.60.bam"
samtools view -h $caps_plus |awk '$0~/^@/||$3~/chr/'|samtools view -bS - > $chrbam
samtools view -s 0.60 -q 10 -b $chrbam > $down
module load Python/3.7.4-GCCcore-8.3.0
/gpfs2/well/ludwig/users/ebu571/Python/asTair_p3.7.4-skylake/bin/astair call -sc True \
        -i $down --known_snp /users/ludwig/ebu571/ebu571/CAPS_pub/resource/SNV.vcf.gz -f /users/ludwig/ebu571/ebu571/CAPS_pub/resource/mm9_nochrM.fa \
        -mq 10 \
        -t 8 -m mCtoT -co CpG \
        --gz -d astair_output/
module load BEDTools/2.30.0-GCC-10.2.0
zcat astair_output/CAPS_mESC_merged_bam.md.q10.chr_0.60_mCtoT_CpG.mods.gz | sort -k1,1 -k2,2n | bedtools intersect -v -a - -b resource/mm9-blacklist.bed > astair_output/CAPS_mESC_merged_bam.md.q10.chr_0.60_mCtoT_CpG.rmblack.bedGraph

cat astair_output/CAPS_mESC_merged_bam.md.q10.chr_0.60_mCtoT_CpG.rmblack.bedGraph | \
                 awk '$11=="No" && $5+$6>0' | sort -k 1,1 -k2,2n | \
                 bedtools intersect -a resource/cpgIslandExt.mm9.4kbflanking.sorted.bed -b - -wa -wb -sorted > meth/CAPS_mESC_merged_bam.md.q10.chr_0.60_mCtoT_CpG.rmblack_CGI.bed


caps="/users/ludwig/ebu571/ebu571/CAPS_pub/bedtool/GSM4708554_caps_mm9_mESC_CpG.rmblacklist.rmsnv.bed"
ace="/users/ludwig/ebu571/ebu571/CAPS_pub/bedtool/GSE116016_ACE-Seq_WT.ESC.mm9_CG.rmblacklist.rmsnv.bed"

cat $caps | awk '$5+$6>0' | sort -k 1,1 -k2,2n | bedtools intersect -a resource/cpgIslandExt.mm9.4kbflanking.sorted.bed -b - -wa -wb -sorted > meth/GSM4708554_caps_mm9_mESC_CpG.rmblacklist.rmsnv_CGI.bed
cat $ace | awk '$5>0' | sort -k 1,1 -k2,2n | bedtools intersect -a resource/cpgIslandExt.mm9.4kbflanking.sorted.bed -b - -wa -wb -sorted > meth/GSE116016_ACE-Seq_WT.ESC.mm9_CG.rmblacklist.rmsnv_CGI.bed