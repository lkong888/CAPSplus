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

# Log locations which are relative to the current
# working directory of the submission
#$ -o tidy_igv/IGV_output_ref.log
#$ -e tidy_igv/IGV_error_ref.log

# bedgraph to bigwig
#cat /users/ludwig/ebu571/ebu571/CAPS_pub/resource/mm9_size.txt | sort -k1,1 -k2,2n >  /users/ludwig/ebu571/ebu571/CAPS_pub/resource/mm9.chrom.sizes

##CAPS plus data
for i in `ls astair_output/*CpG.chr_rmblack.bedGraph`
do
(cat $i | grep ^chr[1-9,X,Y] | awk -F "\t" '$11=="No" && $5+$6>5' | awk '{{OFS="\t"}}{{print $1, $2, $3, $4}}' | sort -k1,1 -k2,2n)  > $i.sort.bed
/well/ludwig/users/ebu571/conda/skylake/envs/env_bedtobw/bin/bedGraphToBigWig  $i.sort.bed /users/ludwig/ebu571/ebu571/CAPS_pub/resource/mm9.chrom.sizes $i.bw 
mv $i.bw tidy_igv/
done

#load CAPS, TAB-seq and ACE-seq data
caps=/users/ludwig/ebu571/ebu571/CAPS_pub/astair_output/CAPS_bam.md.q10_mCtoT_CpG.nochrM.mods.rmblack.bed
tab=/users/ludwig/ebu571/ebu571/CAPS_pub/bedtool/TAB_bam.md.chr.meth.CpG.rmblacklist.rmsnv.bed
ace=/users/ludwig/ebu571/ebu571/CAPS_pub/bedtool/GSE116016_ACE-Seq_WT.ESC.mm9_CG.rmblacklist.rmsnv.bed
tab_brain=/users/ludwig/ebu571/ebu571/CAPS_haiqi/jf_meth/GSE46710_Ad_Front.hmC_sites_FDR_0.01.hg38.bed

cat $caps | grep ^chr[1-9,X,Y] | awk -F "\t" '$11=="No" && $5+$6>5' | awk '{{OFS="\t"}}{{print $1, $2, $3, $4}}' | sort -k1,1 -k2,2n  > tidy_igv/caps_pub.chr.filtered.bed
cat $tab | grep ^chr[1-9,X,Y] | awk -F "\t" '$5+$6>5' |  awk '{{OFS="\t"}}{{print $1, $2, $3, $4/100}}'  | sort -k1,1 -k2,2n  > tidy_igv/Tab.chr.filtered.bed
cat $ace | grep ^chr[1-9,X,Y] | awk -F "\t" '$5>5' | awk '{{OFS="\t"}}{{print $1, $2, $3, $4/$5}}' | sort -k1,1 -k2,2n  > tidy_igv/GSE116016_ACE-Seq.chr.filtered.bed
cat $tab_brain | awk '$8!="Fail"&& $12>=5' | grep ^chr[1-9,X,Y] | awk '{{OFS="\t"}}{{print $9, $10, $11, $13/$12}}' | sort -k1,1 -k2,2n -u > tidy_igv/GSE46710_Ad_Front.hmC_sites_FDR_0.01.hg38.chr.filtered.bed

/well/ludwig/users/ebu571/conda/skylake/envs/env_bedtobw/bin/bedGraphToBigWig tidy_igv/caps_pub.chr.filtered.bed /users/ludwig/ebu571/ebu571/CAPS_pub/resource/mm9.chrom.sizes tidy_igv/caps_pub.chr.filtered.bed.bw
/well/ludwig/users/ebu571/conda/skylake/envs/env_bedtobw/bin/bedGraphToBigWig tidy_igv/Tab.chr.filtered.bed /users/ludwig/ebu571/ebu571/CAPS_pub/resource/mm9.chrom.sizes tidy_igv/Tab.chr.filtered.bed.bw
/well/ludwig/users/ebu571/conda/skylake/envs/env_bedtobw/bin/bedGraphToBigWig tidy_igv/GSE116016_ACE-Seq.chr.filtered.bed /users/ludwig/ebu571/ebu571/CAPS_pub/resource/mm9.chrom.sizes tidy_igv/GSE116016_ACE-Seq.chr.filtered.bed.bw
/well/ludwig/users/ebu571/conda/skylake/envs/env_bedtobw/bin/bedGraphToBigWig tidy_igv/GSE46710_Ad_Front.hmC_sites_FDR_0.01.hg38.chr.filtered.bed /users/ludwig/ebu571/ebu571/CAPS_haiqi/resource/hg38.chrom.sizes tidy_igv/GSE46710_Ad_Front.bed.bw

