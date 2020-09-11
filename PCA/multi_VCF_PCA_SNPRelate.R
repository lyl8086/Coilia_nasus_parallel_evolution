#!/usr/bin/R

library(ggpubr)
library(gtools)
library(SNPRelate)
library(SeqVarTools)

argv    <- commandArgs(TRUE)
vcf_dir <- argv[1]
file_names <- list.files(path=vcf_dir, pattern='.vcf$')
file_names <- mixedsort(file_names)
plots1 <- list()
plots2 <- list()
a=ceiling(length(file_names)/3)

for (i in 1:length(file_names)) {
    vcf <- file_names[i]
	cat("  Processing vcf", vcf, "...")
    LG  <- gsub(".vcf","",vcf)
    seqVCF2GDS(paste(vcf_dir,vcf,sep="/"), ".tmp.gds", parallel=8, verbose=FALSE)
	gds <- seqOpen(".tmp.gds")
	pc  <- snpgdsPCA(gds, eigen.cnt=length(seqGetData(gds, "sample.id")),
		autosome.only=F, remove.monosnp=F, algorithm="randomized", verbose=FALSE)
	pc.va <- pc$varprop*100
    
	# pca dataframe.
	dat <- data.frame(
		"ind"=pc$sample.id, 
		"Populations"=substr(pc$sample.id,1,2),
		"PC1"=pc$eigenvect[,1], "PC2"=pc$eigenvect[,2],
		stringsAsFactors = FALSE)
        
	# pca scatter plot.
    lab1<-paste("PC1"," (",sprintf("%.2f%%", pc.va[1]),")", sep="")
    lab2<-paste("PC2"," (",sprintf("%.2f%%", pc.va[2]),")", sep="")
    p<-ggscatter(dat,x="PC1",y="PC2",size=2, font.label = c(0.1, "plain"),
    color="Populations",shape="Populations",
    title=LG,xlab=lab1,ylab=lab2) +
    theme(
    legend.title=element_blank(),
	legend.position="left",
    plot.title=element_text(hjust=0.5))
    
    if (i%%3 != 1) {p <- p + theme(legend.position="none")}
    plots1[[i]] <- p
    # clean files.
	seqClose(gds)
	cat("done.\n")
}

pdf("Multi_PCA.SNPRelate.pdf", width=14, height=a*3)
ggarrange(plotlist=plots1, ncol = 3, nrow=a) 
dev.off()

