# Modular oxidation of cytosine modifications and their application in direct and quantitative sequencing of 5-hydroxymethylcytosine
Authors: Haiqi Xu,1,2,3 Jinfeng Chen,1,2,* Jingfei Cheng,1,2,* Linzhen Kong,1,2,* Xiufei Chen,1,2 Masato Inoue,1,2 Yibin Liu,4,5 Meiping Zhao,3,† Chun-Xiao Song1,2,

# Introduction
These scripts are for CAPSplus data analysis. The following steps are included:
1. Data preprocessing: fastp, alignment, deduplicate
2. Methylation calling by asTair. Check conversion rate on spike-in.
3. Remove blacklisted regions and known SNV. 
4. Downstream analysis: correlation, genomic anotations, cnv analysis and IGV visualisation.

# Spike in
CpG-methylated lambda DNA, 2 kb unmodified spike-in, 144-mer synthetic 5hmC spike-in

# Reference
mm9 for mESCs samples; hg38 for brain samples

# 1. Data preprocessing
 astair_spikein.smk   ###check conversion rate and false positive
 
 standard_caps_mESC.smk   ###alignment, mark duplicate, calling methylation, filtering for mESCs
 
 standard_caps_brain.smk   ### alignment, mark duplicate, calling methylation, filtering for brain samples

# 2. Raw signals correlation between CAPSplus, CAPS, TAB-seq and ACE-seq(Figure 4c, Figure S13)
caps_filtering.r

raw_signal.sh

figure_mESC.r

# 3. Raw signal correlation between mESCs replicates (Figure S12)
mESC_replicate.r

# 4. Statistical testing of high confidence 5hmC sites and genomic annotations(Figure S16)
 annotation.sh
 
 CAPS_annotation.r
 
 # 5. Subsampling and coverage analysis(Figure S15)
  SI12.sh
  
  figure_mESC.r
  
 # 6. Downstream analysis for normal brain and glioblastoma (Figure 4d,f, Figure S17)
 brain_work.sh
 
 figure_brain.r
 
 # 8. Per-base quality plot(Figure S14)
 phred_score.sh
 
 # 9. IGV visualisation(Figure 4e)
 bedtobw.sh
 
 
 ####Note: the data were co-analysed by Jingfei Cheng and Linzhen Kong
  
  
  
  
  
  
