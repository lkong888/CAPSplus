#####raw signal comparision between healthy and GBM sample.
.libPaths("/well/ludwig/users/ebu571/R/4.0/skylake")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
library(readr)
options(bitmapType='cairo-png')
setwd("/users/ludwig/ebu571/ebu571/CAPS_haiqi/jf_meth")
##load data
meth1 <- fread("brain_healthy_29Jul2022_S3.md.astair_CpG.filtered.bedGraph",  header = FALSE, nThread = 8, stringsAsFactors = FALSE, data.table = FALSE, sep = "\t", verbose = T)
meth2 <- fread("brain_tumor_29Jul2022_S4.md.astair_CpG.filtered.bedGraph",  header = FALSE, nThread = 8, stringsAsFactors = FALSE, data.table = FALSE, sep = "\t", verbose = T)
colnames(meth1) <- c("chr","start","end","meth","mC","aC")
colnames(meth2) <- c("chr","start","end","meth","mC","aC")

chrs <- paste("chr", c(1:19, "X","Y"), sep = "")
depth <- 0
meth1<- meth1[which(meth1$chr %in% chrs),]
meth1 <- meth1[meth1$aC>depth,]
meth2 <- meth2[which(meth2$chr%in%chrs),]
meth2 <- meth2[meth2$aC>depth,]
meth <- merge(meth1[,c(1:6)], meth2[,c(1:6)], by=c("chr","start","end")) 
binsize <- 1000
meth$bin <- floor(meth$start/binsize)
meth <- meth %>% dplyr::select("chr","bin","meth.x", "mC.x", "aC.x", "meth.y", "mC.y", "aC.y")

meth_agg <- aggregate(cbind(mC.x, aC.x, mC.y, aC.y) ~ bin+chr, meth, mean)

meth_agg$meth_ratio.normal <- meth_agg$mC.x / meth_agg$aC.x
meth_agg$meth_ratio.tumour <- meth_agg$mC.y / meth_agg$aC.y

seldat <- meth_agg  %>% dplyr::select("meth_ratio.normal","meth_ratio.tumour")%>% reshape2::melt()
seldat$variable <- gsub("meth_ratio.normal","healthy",seldat$variable)
seldat$variable <- gsub("meth_ratio.tumour","tumour",seldat$variable)
write.table(meth_agg, "/users/ludwig/ebu571/ebu571/CAPS_haiqi/plots/brain_mod_level_in_1kb_bins.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
mean(meth_agg$meth_ratio.tumour) ##0.02478693
mean(meth_agg$meth_ratio.normal) ##[1] 0.1847834

meth_agg <- read_delim("/users/ludwig/ebu571/ebu571/CAPS_haiqi/plots/brain_mod_level_in_1kb_bins.txt",delim="\t",col_names=FALSE)
colnames(meth_agg) <- c("bin", "chr", "mC.x", "aC.x","mC.y", "aC.y","meth_ratio.normal", "meth_ratio.tumour")
p <- ggplot(seldat,aes(x=variable,y=value, fill=variable)) + geom_violin(trim = FALSE, adjust=4) +
  stat_summary(fun = mean, geom = "point",color = "black") +
  stat_summary(fun = mean, geom="text", aes(label = round(..y.., 2)), hjust = -0.2, vjust=-1.5)+
  theme_light() + xlab("") +  theme(axis.text = element_text(color = "black")) +
  ylab("hmCpG modification level in 1kb") +
  scale_fill_manual(values=c("#67A9CF","#EF8A62"))+
  theme(legend.position = "none") 
ggsave("/users/ludwig/ebu571/ebu571/CAPS_haiqi/plots/hmCG_1kb_brain.pdf",p, width=3, height = 3)
################################################################################################


#######################################################################

#meth1 <- fread("GSE46710_Ad_Front.hmC_sites_FDR_0.01.hg38.bed", header = FALSE, nThread = 8, stringsAsFactors = FALSE, data.table = FALSE, sep = "\t", verbose = T, drop=c(1:8), fill=TRUE)
meth1 <- fread("GSE46710_Ad_Front.hmC_sites_FDR_0.01.hg38_filtered.bed", header = FALSE, nThread = 8, stringsAsFactors = FALSE, data.table = FALSE, sep = "\t", verbose = T, fill=TRUE)
meth2 <- fread("brain_healthy_29Jul2022_S3.md.astair_CpG.filtered.chr.bedGraph", header = FALSE, nThread = 8, stringsAsFactors = FALSE, data.table = FALSE, sep = "\t", verbose = T)
colnames(meth1) <- c("chr","start","end","aC","mC","meth")
colnames(meth2) <- c("chr","start","end","meth","mC","aC")
meth <- merge(meth1[,c(1:6)], meth2[,c(1:6)], by=c("chr","start","end"))
binsize <- 10000
meth$bin <- floor(meth$start/binsize)
meth_agg <- aggregate(.~bin+chr, meth[meth$chr %in% paste0("chr",c(seq(1,22),"X","Y")),] %>% 
                        dplyr::select("chr","bin","meth.x", "mC.x", "aC.x", "meth.y", "mC.y", "aC.y"), sum)

write.table(meth_agg, "/users/ludwig/ebu571/ebu571/CAPS_haiqi/plots/brain_pub.10k.cor_1.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
meth_agg <- fread("/users/ludwig/ebu571/ebu571/CAPS_haiqi/plots/brain_pub.10k.cor_1.txt", header = F, nThread = 8, stringsAsFactors = F, data.table = F, sep = "\t", verbose = T)
colnames(meth_agg) <- c("chr","bin","meth.x", "mC.x", "aC.x", "meth.y", "mC.y", "aC.y")
meth_agg$meth_ratio.pub <- meth_agg$mC.x / meth_agg$aC.x
meth_agg$meth_ratio.healthy <- meth_agg$mC.y / meth_agg$aC.y
cov_min <- 500; cov_max <- 5000
meth_agg[meth_agg$aC.x>cov_min & meth_agg$aC.x<cov_max & meth_agg$aC.y>cov_min & meth_agg$aC.y<cov_max,  ] %>% dplyr::select("meth_ratio.pub","meth_ratio.healthy") %>% cor(method="pearson")
# 0.8899921 ###0.890582
meth_agg[meth_agg$aC.x>cov_min & meth_agg$aC.x<cov_max & meth_agg$aC.y>cov_min & meth_agg$aC.y<cov_max,  ] %>% dplyr::select("meth_ratio.pub","meth_ratio.healthy") %>% cor(method="spearman")
# 0.8837846  ###0.8842603
meth_agg_sel <- meth_agg[meth_agg$aC.x>cov_min & meth_agg$aC.x<cov_max & meth_agg$aC.y>cov_min & meth_agg$aC.y<cov_max,  ] 
pdf("/users/ludwig/ebu571/ebu571/CAPS_haiqi/plots/brain_pub.10k.cor.pdf", width = 3.5, height = 4)
smoothScatter(x=meth_agg_sel$meth_ratio.healthy, y= meth_agg_sel$meth_ratio.pub, xlab=("caps\nhealthy_brain( 10k cor:0.89)"), ylab=("tab-seq"), xlim=c(0.00,0.80), ylim=c(0.00,0.80))
dev.off()
###########################################
#test 
meth1 <- fread("merged_bam.sort.deduplicated_CpG_filtered.bed", header = FALSE, nThread = 8, stringsAsFactors = FALSE, data.table = FALSE, sep = "\t", verbose = T, fill=TRUE)
meth2 <- fread("brain_healthy_29Jul2022_S3.md.astair_CpG.filtered.chr.bedGraph", header = FALSE, nThread = 8, stringsAsFactors = FALSE, data.table = FALSE, sep = "\t", verbose = T)
colnames(meth1) <- c("chr","start","end","meth","mC","uC")
colnames(meth2) <- c("chr","start","end","meth","mC","aC")
meth11 <- meth1
meth1 <- meth11
meth1 <- meth1[(meth1$mC+meth1$uC)>=5,]
meth <- merge(meth1[,c(1:6)], meth2[,c(1:6)], by=c("chr","start","end"))
binsize <- 10000
meth$bin <- floor(meth$start/binsize)
meth_agg <- aggregate(.~bin+chr, meth[meth$chr %in% paste0("chr",c(seq(1,22),"X","Y")),] %>% 
                        dplyr::select("chr","bin","meth.x", "mC.x", "uC", "meth.y", "mC.y", "aC"), sum)

meth_agg$meth_ratio.pub <- meth_agg$mC.x / (meth_agg$mC.x + meth_agg$uC)
meth_agg$meth_ratio.healthy <- meth_agg$mC.y / meth_agg$aC
cov_min <- 500; cov_max <- 5000
meth_agg[(meth_agg$mC.x + meth_agg$uC)>cov_min & (meth_agg$mC.x + meth_agg$uC)<cov_max & meth_agg$aC>cov_min & meth_agg$aC<cov_max,  ] %>% dplyr::select("meth_ratio.pub","meth_ratio.healthy") %>% cor(method="pearson")
# 0.8899921 ###0.890582 ###0.9139511

write.table(meth_agg, "/users/ludwig/ebu571/ebu571/CAPS_haiqi/plots/brain_pub.10k.cor_2.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
meth_agg[(meth_agg$mC.x + meth_agg$uC)>cov_min & (meth_agg$mC.x + meth_agg$uC)<cov_max & meth_agg$aC>cov_min & meth_agg$aC<cov_max,  ] %>% dplyr::select("meth_ratio.pub","meth_ratio.healthy") %>% cor(method="spearman")
# 0.8837846  ###0.8842603
meth_agg_sel <- meth_agg[(meth_agg$mC.x + meth_agg$uC)>cov_min &(meth_agg$mC.x + meth_agg$uC)<cov_max & meth_agg$aC>cov_min & meth_agg$aC<cov_max,  ] 
pdf("/users/ludwig/ebu571/ebu571/CAPS_haiqi/plots/brain_pub.10k.cor.pdf", width = 3.5, height = 4)
smoothScatter(x=meth_agg_sel$meth_ratio.healthy, y= meth_agg_sel$meth_ratio.pub, xlab=("caps\nhealthy_brain( 10k cor:0.91)"), ylab=("tab-seq"), xlim=c(0.00,0.80), ylim=c(0.00,0.80))
dev.off()
