setwd("/gpfs2/well/ludwig/users/ebu571/CAPS_haiqi")
.libPaths("/well/ludwig/users/ebu571/R/4.0/skylake")
library(data.table)
library(dplyr)
library(Smisc)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(grid)
library(parallel)
library(MASS)

chrs <- paste("chr", c(1:19, "X","Y"), sep = "")
cutoff <- 0

##CpG sites
cpg1 <- read_delim("astair_output/mESC_1_29Jul2022_S1.md.q10_mCtoT_CpG.chr_rmblack.bedGraph.gz", delim="\t", col_names=FALSE)
colnames(cpg1) <- c("chr","start","end","mod_level","mod","unmod", "ref", "alt", "specific_context","context","snv","total_depth")
nrow(cpg1)##[1] 40509751
cpg1 <- cpg1[which(cpg1$chr %in% chrs),]
cpg1<- cpg1[cpg1$mod+cpg1$unmod>=cutoff & cpg1$snv=="No",]  %>% dplyr::select(chr, start, end, mod, unmod)
nrow(cpg1)###[1] 40064449

cpg2 <- read_delim("astair_output/mESC_2_29Jul2022_S2.md.q10_mCtoT_CpG.chr_rmblack.bedGraph.gz", delim="\t", col_names=FALSE)
colnames(cpg2) <- c("chr","start","end","mod_level","mod","unmod", "ref", "alt", "specific_context","context","snv","total_depth")
#40528646
cpg2 <- cpg2[which(cpg2$chr %in% chrs),]
cpg2 <- cpg2[cpg2$mod+cpg2$unmod>=cutoff & cpg2$snv=="No",] %>% dplyr::select(chr, start, end, mod, unmod)
nrow(cpg2)##[1] 40083325

dat1 <- merge(cpg1, cpg2, by=c("chr", "start", "end"), all.x=T, all.y=T)
colnames(dat1) <- c("chr", "start", "end", "r1_mod", "r1_unmod", "r2_mod", "r2_unmod")

table(!is.na(dat1$r1_mod))
#FALSE     TRUE 
#203823 40064449  
table(!is.na(dat1$r2_mod))
#FALSE     TRUE 
#184947 40083325
table(!is.na(dat1$r1_mod) & !is.na(dat1$r2_mod))
#FALSE     TRUE 
#388770 39879502 

library(VennDiagram)
pdf("plot/astair_mESC_replicate_venn.pdf")
draw.pairwise.venn(area1=40064449, area2=40083325,cross.area=39879502,
                   category=c("mESC1","mESC2"),fill=c("aquamarine","grey"))
dev.off()




###read data: 
##1. rawa signal
caps_r1 <- read_delim("astair_output/mESC_1_29Jul2022_S1.md.q10_mCtoT_CpG.chr_rmblack_raw_signal.bed", delim="\t", col_names=FALSE)
colnames(caps_r1) <- c("chr", "start", "end", "mod", "unmod")
caps_r1 <- caps_r1[which(caps_r1$chr %in% chrs),]


caps_r2 <- read_delim("astair_output/mESC_2_29Jul2022_S2.md.q10_mCtoT_CpG.chr_rmblack_raw_signal.bed", delim="\t", col_names=FALSE)
colnames(caps_r2) <- c("chr", "start", "end", "mod", "unmod")
caps_r2 <- caps_r2[which(caps_r2$chr %in% chrs),]

rep <- merge(caps_r1, caps_r2, by=c("chr", "start", "end"))
colnames(rep) <- c("chr", "start", "end", "r1_mod", "r1_unmod", "r2_mod", "r2_unmod")

selrep <- rep[!(is.na(rep$r1_mod) & is.na(rep$r1_unmod)&
                is.na(rep$r2_mod) & is.na(rep$r2_unmod)),]
nrow(selrep)
#[1] 254416

selrep$r1_signal <- selrep$r1_mod/(selrep$r1_mod+selrep$r1_unmod)
selrep$r2_signal <- selrep$r2_mod/(selrep$r2_mod+selrep$r2_unmod)
options(bitmapType='cairo-png')
hist((selrep$r1_mod+selrep$r1_unmod), breaks=1000, xlim = c(0, 6000))
hist((selrep$r2_mod+selrep$r2_unmod), breaks=1000, xlim = c(0, 6000))

depth_cutoff <- 500
selrep2 <- selrep[(selrep$r1_mod+selrep$r1_unmod)>=depth_cutoff & (selrep$r1_mod+selrep$r1_unmod) <= 5000 &
                    (selrep$r2_mod+selrep$r2_unmod)>=depth_cutoff & (selrep$r2_mod+selrep$r2_unmod) <= 5000,]
nrow(selrep2)#[1] 206778


###reproduction
idx_rep <- !is.na(selrep2$r1_signal) & !is.na(selrep2$r2_signal)
pdf("plot/astair_replicate1_vs_2_raw_signal.pdf")
smoothScatter(selrep2$r1_signal[idx_rep], selrep2$r2_signal[idx_rep], xlab = "mESC1", ylab = "mESC2", xlim = c(0, 0.25), ylim = c(0,0.25))
dev.off()
cor(selrep2$r1_signal[idx_rep], selrep2$r2_signal[idx_rep]) ###, use = "complete.obs"
##[1] 0.8505002
cor(selrep2$r1_signal[idx_rep], selrep2$r2_signal[idx_rep], method = "spearman") ###, use = "complete.obs"
##[1] 0.8180385

