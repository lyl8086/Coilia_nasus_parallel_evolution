library(ggpubr)
library(gtools)
library(reshape2)
library(SeqVarTools)

argv    <- commandArgs(TRUE)
vcf_dir <- argv[1]
file_names <- list.files(path=vcf_dir, pattern='.vcf$')
file_names <- mixedsort(file_names)
plots1 <- list()
plots2 <- list()
a <- length(file_names)

pca_scatter <- function(pc=pc, va=pc.va, x="PC1", y="PC2", v1=1, v2=2){
	
	# create pca dataframe.
	dat <- data.frame(
		"ind"=row.names(pc$eigenvect), 
		"Populations"=substr(row.names(pc$eigenvect),1,2),
		x=pc$eigenvect[,v1], y=pc$eigenvect[,v2],
		stringsAsFactors = FALSE)
	colnames(dat)[3:4] <- c(x, y)
	lab1<-paste(x," (",sprintf("%.2f%%", 100*va[v1]/sum(va)),")", sep="")
    lab2<-paste(y," (",sprintf("%.2f%%", 100*va[v2]/sum(va)),")", sep="")
    p <- ggscatter(dat,x=x,y=y,size=2, font.label = c(0.1, "plain"),
		color="Populations",shape="Populations",
		title=LG,xlab=lab1,ylab=lab2) +
		theme(
		legend.title=element_blank(),
		legend.position="left",
		plot.title=element_text(hjust=0.5))
	return(p)
}

for (i in 1:length(file_names)) {
    vcf <- file_names[i]
	cat("  Processing vcf", vcf, "...")
    LG  <- gsub(".vcf","",vcf)
    seqVCF2GDS(paste(vcf_dir,vcf,sep="/"), ".tmp.gds", parallel=8, verbose=FALSE)
	gds <- seqOpen(".tmp.gds")
	pc  <- pca(gds, eigen.cnt=length(seqGetData(gds, "sample.id")))
	pc.va <- pc$eigenval
    plots1[[i]] <- ggarrange(
		pca_scatter(pc, pc.va, "PC1", "PC2", 1, 2),
		pca_scatter(pc, pc.va, "PC1", "PC3", 1, 3),
		pca_scatter(pc, pc.va, "PC2", "PC3", 2, 3),
		nrow = 1	
	)
    # clean files.
	seqClose(gds)
	cat("done.\n")
}

pdf("Multi_PCA.SeqVarTools.pdf", width=14, height=a*3)
ggarrange(plotlist=plots1, ncol = 1)
dev.off()

