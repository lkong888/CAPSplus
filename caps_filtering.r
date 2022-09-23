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
depth <- 0
#current data:
caps_merge <- read_delim("astair_output/CAPS_mESC_merged_bam.md.q10_mCtoT_CpG.chr_rmblack.bedGraph.gz", delim="\t", col_names=FALSE)
colnames(caps_merge) <- c("chr","start","end","mod_level","mod","unmod", "ref", "alt", "specific_context","context","snv","total_depth")
caps_merge <- caps_merge[caps_merge$snv == "No",]
caps_merge <- caps_merge[,c("chr","start","end","mod_level","mod","unmod","context","total_depth")]
nrow(caps_merge)
caps_merge<- caps_merge[which(caps_merge$chr %in% chrs),] ###rm chrM
nrow(caps_merge)
caps_merge <- caps_merge[(caps_merge$mod+caps_merge$unmod)>depth,]
nrow(caps_merge)

###caps_yibin liu
caps_pub <- read_delim("/users/ludwig/ebu571/ebu571/CAPS_pub/bedtool/GSM4708554_caps_mm9_mESC_CpG.rmblacklist.rmsnv.bed", delim="\t", col_names=FALSE)
colnames(caps_pub) <- c("chr", "start", "end", "mod_level", "mod", "unmod")
nrow(caps_pub)
caps_pub<- caps_pub[which(caps_pub$chr %in% chrs),]
nrow(caps_pub)
caps_pub <- caps_pub[caps_pub$mod+caps_pub$unmod>depth,]
nrow(caps_pub)

###published ace_seq 
ace_pub <- read_delim("/users/ludwig/ebu571/ebu571/CAPS_pub/bedtool/GSE116016_ACE-Seq_WT.ESC.mm9_CG.rmblacklist.rmsnv.bed", delim="\t", col_names=FALSE)
colnames(ace_pub) <- c("chr", "start", "end", "mod", "aC", "strand")
nrow(ace_pub)
ace_pub <- ace_pub[which(ace_pub$chr %in% chrs),]
nrow(ace_pub)
ace_pub <- ace_pub[ace_pub$aC>depth,]
nrow(ace_pub)

###TAB-seq
tab <- read_delim("/users/ludwig/ebu571/ebu571/CAPS_pub/bedtool/TAB_bam.md.chr.meth.CpG.rmblacklist.rmsnv.bed",  delim="\t", col_names=FALSE)
colnames(tab) <- c("chr", "start", "end", "mod_level", "mod", "unmod")
nrow(tab)
tab <- tab[which(tab$chr %in% chrs),]
nrow(tab)
tab <- tab[tab$mod+tab$unmod>depth,]
nrow(tab)

###merge data

cpg1 <- merge(caps_merge, caps_pub, by=c("chr", "start", "end"))
nrow(cpg1)
colnames(cpg1) <- c("chr", "start", "end", "caps_merge_mod_level", "caps_mereg_mod", "caps_merge_unmod", "caps_merge_context", "caps_merge_total_depth", "caps_pub_mod_level",  "caps_pub_mod", "caps_pub_unmod")
cpg1 <- cpg1 %>% dplyr::select(chr, start, end, caps_mereg_mod, caps_merge_unmod, caps_pub_mod, caps_pub_unmod)

cpg2 <- merge(cpg1, ace_pub, by=c("chr", "start", "end"))
colnames(cpg2) <- c("chr", "start", "end", "caps_mereg_mod", "caps_merge_unmod", "caps_pub_mod", "caps_pub_unmod", "ace_pub_mod", "ace_pub_aC", "ace_pub_strand")
nrow(cpg2)
cpg2 <- cpg2 %>% dplyr::select(chr, start, end, caps_mereg_mod, caps_merge_unmod, caps_pub_mod, caps_pub_unmod, ace_pub_mod, ace_pub_aC)

cpg_all <- merge(cpg2, tab, by=c("chr", "start", "end"))
colnames(cpg_all) <-  c("chr", "start", "end", "caps_mereg_mod", "caps_merge_unmod",  "caps_pub_mod", "caps_pub_unmod", "ace_pub_mod", "ace_pub_aC", "tab_mod_level", "tab_mod", "tab_unmod")
nrow(cpg_all)
cpg_all <- cpg_all  %>% dplyr::select(chr, start, end, caps_mereg_mod, caps_merge_unmod, caps_pub_mod, caps_pub_unmod, ace_pub_mod, ace_pub_aC, tab_mod, tab_unmod)

cpg_all <- format(cpg_all, scientific = FALSE)
temp <- cpg_all$ace_pub_mod %>% as.numeric()
cpg_all$ace_pub_mod<- round(temp)

fwrite(cpg_all, "astair_output/caps_all.bedGraph", sep = "\t", col.names = F, row.names = F, quote = FALSE)
