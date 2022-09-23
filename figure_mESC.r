##module load R/4.0.3-foss-2020b
.libPaths("/well/ludwig/users/ebu571/R/4.0/skylake")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
library(readr)
options(bitmapType='cairo-png')
setwd("/users/ludwig/ebu571/ebu571/CAPS_haiqi")

#####correlation########
chrs <- paste("chr", c(1:19, "X","Y"), sep = "")
###read data:
caps_all_raw <- read_delim("astair_output/caps_all_rawsignals.bed", delim="\t", col_names=FALSE)
colnames(caps_all_raw) <- c("chr", "start", "end", "caps_merge_mod", "caps_merge_unmod", "caps_pub_mod", "caps_pub_unmod", "ace_pub_mod", "ace_pub_aC", "tab_mod", "tab_unmod")
caps_all_raw <- caps_all_raw[which(caps_all_raw$chr %in% chrs),]
nrow(caps_all_raw)##[1] 265500

###calculate raw signal
caps_all_raw$caps_merge_signal <- caps_all_raw$caps_merge_mod/(caps_all_raw$caps_merge_mod+caps_all_raw$caps_merge_unmod)
caps_all_raw$caps_pub_signal <- caps_all_raw$caps_pub_mod/(caps_all_raw$caps_pub_mod+caps_all_raw$caps_pub_unmod)
caps_all_raw$ace_signal <- caps_all_raw$ace_pub_mod/caps_all_raw$ace_pub_aC
caps_all_raw$tab_signal <- caps_all_raw$tab_mod/(caps_all_raw$tab_mod+caps_all_raw$tab_unmod)

options(bitmapType='cairo-png')
hist(caps_all_raw$caps_pub_mod+caps_all_raw$caps_pub_unmod, breaks=4000, xlim = c(0, 10000))
hist(caps_all_raw$caps_merge_aC, breaks=4000, xlim = c(0, 10000))
hist(caps_all_raw$ace_pub_aC, breaks=4000, xlim = c(0, 10000))
hist(caps_all_raw$tab_mod+caps_all_raw$tab_unmod, breaks=4000, xlim = c(0, 10000))

depth_cutoff <- 500
caps_all_raw_2 <- caps_all_raw[(caps_all_raw$caps_merge_mod+caps_all_raw$caps_merge_unmod)>=depth_cutoff & (caps_all_raw$caps_merge_mod+caps_all_raw$caps_merge_unmod) <= 5000 &
                                 (caps_all_raw$caps_pub_mod+caps_all_raw$caps_pub_unmod)>=depth_cutoff & (caps_all_raw$caps_pub_mod+caps_all_raw$caps_pub_unmod) <= 5000 &
                                 caps_all_raw$ace_pub_aC >= depth_cutoff & caps_all_raw$ace_pub_aC <= 5000 &
                                 (caps_all_raw$tab_mod+caps_all_raw$tab_unmod) >= depth_cutoff & (caps_all_raw$tab_mod+caps_all_raw$tab_unmod) <= 5000, ]
nrow(caps_all_raw_2)

###scatter plot and correlation
###caps_vs_caps_pub
idx_caps <- !is.na(caps_all_raw_2$caps_merge_signal) & !is.na(caps_all_raw_2$caps_pub_signal)
pdf("plot/astair_caps_vs_caps_pub_raw_signal.pdf")
smoothScatter(caps_all_raw_2$caps_merge_signal[idx_caps], caps_all_raw_2$caps_pub_signal[idx_caps], xlab = "CAPS+ raw signal", ylab = "CAPS raw signal", xlim = c(0, 0.2), ylim = c(0,0.2))
dev.off()
cor(caps_all_raw_2$caps_merge_signal[idx_caps],  caps_all_raw_2$caps_pub_signal[idx_caps]) 
##[1] 0.8701933 astair

cor(caps_all_raw_2$caps_merge_signal[idx_caps],  caps_all_raw_2$caps_pub_signal[idx_caps], method = "spearman")
##[1] 0.852222 astair


###caps_merge_vs ace_pub
idx_ace <- !is.na(caps_all_raw_2$caps_merge_signal) & !is.na(caps_all_raw_2$ace_signal)
pdf("plot/astair_caps_vs_ace_pub_raw_signal.pdf")
smoothScatter(caps_all_raw_2$caps_merge_signal[idx_ace], caps_all_raw_2$ace_signal[idx_ace], xlab="CAPS+ raw signal", ylab="ACE-seq raw signal", xlim = c(0, 0.2), ylim = c(0,0.2))
dev.off()
cor(caps_all_raw_2$caps_merge_signal[idx_ace], caps_all_raw_2$ace_signal[idx_ace]) ##[1] 0.6362912 astair
cor(caps_all_raw_2$caps_merge_signal[idx_ace], caps_all_raw_2$ace_signal[idx_ace], method = "spearman") ##[1] 0.6883624 astair


###caps_merge_vs_tab
idx_tab <- !is.na(caps_all_raw_2$caps_merge_signal) & !is.na(caps_all_raw_2$tab_signal)
pdf("plot/astair_caps_merge_vs_tab_raw_signal.pdf")
smoothScatter(caps_all_raw_2$caps_merge_signal[idx_tab], caps_all_raw_2$tab_signal[idx_tab], xlab="CAPS+ raw signal", ylab="TAB-seq raw signal", xlim = c(0, 0.2), ylim = c(0,0.2))
dev.off()
cor(caps_all_raw_2$caps_merge_signal[idx_tab], caps_all_raw_2$tab_signal[idx_tab])##[1] 0.7928195 astair
cor(caps_all_raw_2$caps_merge_signal[idx_tab], caps_all_raw_2$tab_signal[idx_tab], method = "spearman") ##[1] 0.8048322 astair


###TAB_vs_ace
idx_n <- !is.na(caps_all_raw_2$tab_signal) & !is.na(caps_all_raw_2$ace_signal)
pdf("plot/astair_tab_vs_ace_raw_signal.pdf")
smoothScatter(caps_all_raw_2$tab_signal[idx_n], caps_all_raw_2$ace_signal[idx_n], xlab="TAB-seq raw signal", ylab="ACE-seq raw signal", xlim = c(0, 0.2), ylim = c(0,0.2))
dev.off()
cor(caps_all_raw_2$tab_signal[idx_n], caps_all_raw_2$ace_signal[idx_n]) ##[1] 0.7596591 astair
cor(caps_all_raw_2$tab_signal[idx_n], caps_all_raw_2$ace_signal[idx_n], method = "spearman")###[1] 0.7864829 astair
############################## correlation analysis between CAPSplus and CAPS, ACE-seq, TAB-seq 


#######coverage over CGI#####
capsplus <- read_delim("meth/CAPS_mESC_merged_bam.md.q10.chr_0.60_mCtoT_CpG.rmblack_CGI.bed",delim="\t",col_names=FALSE)
colnames(capsplus) <- c("gchr", "gstart", "gend", "index", "chr","start","end","mod_level","mod","unmod", "ref", "alt", "specific_context","context","snv","total_depth")
capsplus$aC<- capsplus$mod+capsplus$unmod
capsplus <- capsplus[,c("index", "aC")]

caps<- read_delim("meth/GSM4708554_caps_mm9_mESC_CpG.rmblacklist.rmsnv_CGI.bed",delim="\t",col_names=FALSE)
colnames(caps) <- c("gchr", "gstart", "gend", "index","chr", "start", "end", "mod_level", "mod", "unmod")
caps$aC <- caps$mod+caps$unmod
caps <- caps[,c("index", "aC")]

ace<- read_delim("meth/GSE116016_ACE-Seq_WT.ESC.mm9_CG.rmblacklist.rmsnv_CGI.bed",delim="\t",col_names=FALSE)
colnames(ace) <- c("gchr", "gstart", "gend", "index","chr", "start", "end", "mod", "aC", "strand")
ace <- ace[,c("index", "aC")]

sel_r <- aggregate(aC ~ index,capsplus,mean)
sel_caps <- aggregate(aC ~ index,caps,mean)
sel_ace <- aggregate(aC ~ index, ace, mean)

all_CpG <- merge(sel_r, sel_caps, by.x=c("index"), by.y=c("index"))
colnames(all_CpG) <- c("index", "r_aC", "caps_aC")
all_CpG <- merge(all_CpG, sel_ace, by.x=c("index"), by.y=c("index"))
colnames(all_CpG) <- c("index", "r_aC", "caps_aC", "ace_aC")

mel_CpG <- melt(all_CpG, id.vars =c("index"))

p1 <- ggplot(mel_CpG, aes(x=index, y=value, color=variable)) +
  geom_line()+
  theme(legend.position = "bottom") +
  ylim(5.5,10.5) +
  ylab("coverage") +
  xlab("around CGI") +
  scale_x_continuous(breaks=c(0,20,30,50),
                     labels=c("CpG-4k", "                CpG island","", "CpG+4k"))+
  scale_color_manual(values=c("#CC6600","#446db4","#808080"),name="")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
ggsave("plots/CpG coverage around CGI_0.60.pdf", p1, width=5,height = 3)

