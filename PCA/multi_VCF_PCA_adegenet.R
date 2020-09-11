#!/usr/bin/R

library(ggpubr)
library(gtools)
library(vcfR)
library(adegenet)

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
    LG  <- gsub(".recode.vcf","",vcf)
	vcfr <- read.vcfR(paste(vcf_dir,vcf,sep="/"), verbose=F)
	gldf <- vcfR2genlight(vcfr)
	pop(gldf) <- substr(gldf$ind.names,1,2)
	pca <- glPca(gldf, nf=2, loadings=F, useC=FALSE)
	pc  <- pca$scores
	pc.va <- pca$eig
	pc1 <- pc[,1]
	
	# pca dataframe.
	dat <- data.frame(
		"ind"=gldf$ind.names,
		"Populations"=pop(gldf),
		"PC1"=pc[,1], "PC2"=pc[,2],
		stringsAsFactors = FALSE)
	# pca scatter plot.
    lab1<-paste("PC1"," (",sprintf("%.2f%%", 100*pc.va[1]/sum(pc.va)),")", sep="")
    lab2<-paste("PC2"," (",sprintf("%.2f%%", 100*pc.va[2]/sum(pc.va)),")", sep="")
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
	cat("done.\n")
}

pdf("Multi_PCA.adegenet.pdf", width=14, height=a*3)
ggarrange(plotlist=plots1, ncol = 3, nrow=a) 
dev.off()

