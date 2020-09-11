#!/usr/bin/R

# plink --vcf $in --out out --make-bed --allow-extra-chr --keep-allele-order
# awk -v p=$pop '{n=substr($2,1,p);print $1,$2,$1,n}' out.fam >fam.id
# plink --bfile out --out out --make-grm-bin --keep-allele-order --allow-extra-chr --update-ids fam.id
# gcta64 --grm out --pca --out pca
library(ggpubr)
library(gtools)
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
	
    # pca dataframe.
	dat <- data.frame(
		"ind"=pc[,1],
		"Populations"=pc[,2],
		"PC1"=pc[,3], "PC2"=pc[,4],
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
	system("rm out.* fam.id pca.*")
	cat("done.\n")
}

pdf("Multi_PCA.gcta.pdf", width=14, height=a*3)
ggarrange(plotlist=plots1, ncol = 3, nrow=a) 
dev.off()

