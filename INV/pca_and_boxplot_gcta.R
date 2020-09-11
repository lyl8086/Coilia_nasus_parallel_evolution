#!/usr/bin/R

library(ggpubr)
library(gtools)
library(dplyr)
library(reshape2)
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
    LG  <- gsub(".recode.vcf","",vcf)
	seqVCF2GDS(paste(vcf_dir,vcf,sep="/"), ".tmp.gds", parallel=8, verbose=FALSE)
	gds <- seqOpen(".tmp.gds")
	# vcf to bed.
    system(
        paste("plink --vcf", paste(vcf_dir,vcf,sep="/"),
        "--out out", 
        "--make-bed","--allow-extra-chr",
        "--keep-allele-order"),
		ignore.stdout=T
    )
    system("awk -v p=2 '{n=substr($2,1,p);print $1,$2,$1,n}' out.fam >fam.id")
    # update pop id.
    system(
        paste("plink --bfile out","--out out", 
        "--make-grm-bin","--allow-extra-chr",
		"--mind 0.5",
        "--keep-allele-order","--update-ids fam.id"),,
		ignore.stdout=T
    )
    # pca
    system("gcta64 --grm out --pca --out pca",
		ignore.stdout=T)
    pc <- read.table("pca.eigenvec",header=F, stringsAsFactors = F)
    pc.va <- read.table("pca.eigenval",header=F, stringsAsFactors = FALSE)
	pc.va <- pc.va$V1
	pc1 <- pc[,3]
	names(pc1) <- pc[,1] 
	set.seed(222) # for reproducible
	# partinate to 3 clusters by kmeans.
	clst <- kmeans(pc1, 3, iter.max=100, nstart=100)
	clst.o  <- order(clst$centers)
	hap1 <- names(clst$cluster[which(clst$cluster==clst.o[1])])
	hap2 <- names(clst$cluster[which(clst$cluster==clst.o[2])])
	hap3 <- names(clst$cluster[which(clst$cluster==clst.o[3])])
    write.table(hap1, file=paste0(vcf,".hap1"), quote=F, row.names=F, col.names=F)
	write.table(hap2, file=paste0(vcf,".hap2"), quote=F, row.names=F, col.names=F)
	write.table(hap3, file=paste0(vcf,".hap3"), quote=F, row.names=F, col.names=F)
	# calc het.
	het <- list(
    "hap1"=heterozygosity(gds, margin="by.sample", use.names=T)[hap1], 
	"hap2"=heterozygosity(gds, margin="by.sample", use.names=T)[hap2], 
	"hap3"=heterozygosity(gds, margin="by.sample", use.names=T)[hap3])	
	md <- melt(het, varnames=c("het", "hap"))
	names(md) <- c("het","hap")
	# boxplot.
	hapbox <- ggboxplot(md, x="hap", y="het", 
	color = "hap",add = "jitter", shape="hap", 
	xlab="", title=LG) +
	theme(legend.position="none", plot.title=element_text(hjust=0.5))
	# pca dataframe.
	v <- rep(c("1","2","3"), c(length(hap1),length(hap2),length(hap3)))
	names(v) <- c(hap1,hap2,hap3)
	dat <- data.frame(
		"ind"=pc[,1], 
		"clusters"=v[pc[,1]],
		"Populations"=pc[,2],
		"PC1"=pc[,3], "PC2"=pc[,4],
		stringsAsFactors = FALSE)
	# pca scatter plot.
    lab1<-paste("PC1"," (",sprintf("%.2f%%", 100*pc.va[1]/sum(pc.va)),")", sep="")
    lab2<-paste("PC2"," (",sprintf("%.2f%%", 100*pc.va[2]/sum(pc.va)),")", sep="")
    p<-ggscatter(dat,x="PC1",y="PC2",size=2, font.label = c(0.1, "plain"),
    color="Populations",shape="clusters",
    title=LG,xlab=lab1,ylab=lab2) +
    theme(
    legend.title=element_blank(),
	legend.position="left",
    plot.title=element_text(hjust=0.5))
    
    if (i%%3 != 1) {p <- p + theme(legend.position="none")}
    plots1[[i]] <- p
	plots2[[i]] <- hapbox
    # clean files.
	seqClose(gds)
	system("rm out.* pca.*")
	cat("done.\n")
}

pdf("LDna.PCA.gcta.pdf", width=14, height=a*3)
ggarrange(plotlist=plots1, ncol = 3, nrow=a) 
dev.off()

pdf("LDna.het.gcta.pdf", width=14, height=a*3)
ggarrange(plotlist=plots2, ncol = 3, nrow=a) 
dev.off()
