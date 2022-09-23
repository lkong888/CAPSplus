#module load R/4.1.2-foss-2021b
.libPaths("/well/ludwig/users/ebu571/R/4.1/skylake")
library(CNAclinic)
library(QDNAseq)
library(Biobase)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(pheatmap)
library(ggpubr)
library(RColorBrewer)
library(R.cache)

setwd("/users/ludwig/ebu571/ebu571/CAPS_haiqi/cnv")

bamfiles <- c("brain_healthy_29Jul2022_S3.chr.md.bam", 
              "brain_tumor_29Jul2022_S4.chr.md.bam")
bamnames <- c("healthy", "tumor")
plotcoverage <- function(bamfile, binSize){
    load(paste0(binSize,"hg38.rData"))
    processedData <- processForSegmentation(bamfiles = bamfile, 
        binSize=binSize, 
        userMadeBins = bins, 
        pairedEnds = TRUE,
        chromosomesFilter=c("chrX", "chrY"),           
        cache=FALSE,
        isPaired=TRUE,
        saveCountData=FALSE)

    save(processedData, file=paste0(gsub(".bam","",bamfile),".bin",binSize,".processed.rData"))
    CNAData <- runSegmentation(processedData, 
       genome="hg38",
       segmentType=c("CBS"),
       summaryMethod="mean")
    save(CNAData, file=paste0(gsub(".bam","",bamfile),".bin",binSize,".cna.rData"))

    genomewide_plot <- plotSampleData(CNAData, 
                            showCalls=TRUE, 
                            segmentType="summary",
                            xaxSize=7,
                            mainSize=12)


    genomewide_plot <- genomewide_plot[[1]]
    ggsave(paste0(gsub(".bam",".bin",bamfile),binSize,".genomeplot.pdf"),genomewide_plot)
    #ggsave(paste0(gsub(".bam",".bin",bamfile),binSize,".genomeplot.png"),genomewide_plot)

    chromosomes_plot <- plotSampleData(CNAData, 
        chromosomesFilter=c("X","Y"),
        showCalls=FALSE, 
        segmentType=c("CBS", "PLS", "HMM"),
        pointShape=19,
        pointSize=0.9,   ###0.9
        xaxSize=12,
        pointColor="grey20",
        segHMMColor="hotpink", 
        segCBSColor="skyblue",
        segPLSColor="orange",
        segmentLineWidth=1,
        main="",
        ylim=c(-2, 1.5))

    # A list is returned, accessing 1st element for our single sample
    chromosomes_plot <- chromosomes_plot[[1]]
    ggsave(paste0(gsub(".bam",".bin",bamfile),binSize,".chrplot.pdf"),chromosomes_plot)
    #ggsave(paste0(gsub(".bam",".bin",bamfile),binSize,".chrplot.png"),chromosomes_plot)
}


 for(bamfile in bamfiles){
         plotcoverage(bamfile,500) ###50
 }


 ######test
 #binSize <- 50
 #bins <- createBins(bsgenome=BSgenome.Hsapiens.UCSC.hg38, binSize=binSize)
 #bins$mappability <- calculateMappability(bins,bigWigFile="/users/ludwig/ebu571/ebu571/CAPS_haiqi/cnv/map_hg38_50kb.bw", "/users/ludwig/ebu571/ebu571/bigWigAverageOverBed")
 #wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/bigWigAverageOverBed
 #chmod +x /users/ludwig/ebu571/ebu571/bigWigAverageOverBed
 #bins$blacklist <- calculateBlacklist(bins, bedFiles=c("hg38_blacklist.bed", "hg38_centromere.bed"))

load("brain_tumor_29Jul2022_S4.chr.md.bin500.cna.rData")
genomewide_plot <- plotSampleData(CNAData, 
                            showCalls=TRUE, 
                            segmentType="summary",
                            lossColor="blue",
                            gainColor="red",
                            pointColor="black",
                            segSummaryColor="black",
                            pointShape=19,
                            pointSize=1,
                            xaxSize=7,
                            mainSize=12,
                            segmentLineWidth=1,
                            ylim=c(-2, 2))
genomewide_plot <- genomewide_plot[[1]]
ggsave("brain_tumor_29Jul2022_S4.chr.md.bin500.genomeplot_final.pdf",genomewide_plot, width = 24, height = 4)


load("brain_healthy_29Jul2022_S3.chr.md.bin500.cna.rData")
genomewide_plot <- plotSampleData(CNAData, 
                            showCalls=TRUE, 
                            segmentType="summary",
                            lossColor="blue",
                            gainColor="red",
                            pointColor="grey20",
                            segSummaryColor="black",
                            pointShape=19,
                            pointSize=1,
                            xaxSize=7,
                            mainSize=12,
                            segmentLineWidth=1,
                            ylim=c(-2, 2))
genomewide_plot <- genomewide_plot[[1]]
ggsave("brain_healthy_29Jul2022_S3.chr.md.bin500.genomeplot_final.pdf",genomewide_plot, width = 24, height = 4)
