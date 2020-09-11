argv   <- commandArgs(TRUE)
infile <- argv[1]
out    <- argv[2]

library(LDna)
library(stringi)
library(dplyr)
library(data.table)
# CHR_A         BP_A SNP_A  CHR_B         BP_B SNP_B           R2

cat(paste("Read file", infile, "..."))
x<-fread(infile, header = T)
x<-mutate(x, id1=paste(CHR_A,BP_A,sep="_"),id2=paste(CHR_B,BP_B,sep="_"))
cat("done.\n")

# init matrix
cat("Initiating matrix...")
mat.names<-unique(c(x$id1,tail(x$id2,1)))
n.loci<-length(mat.names)
mat<-matrix(nrow=n.loci, ncol=n.loci, dimnames=list(mat.names,mat.names))
# fill lower tri 
mat[lower.tri(mat)]<-x$R2
cat("done.\n")

# LDna
ldna <- LDnaRaw(mat, digits=3, mc.cores=8)

# get clusters
pdf(paste(out,'.pdf', sep=''), width=7, height=5)
par(mfcol=c(1,2))
clusters0 <- extractClusters(ldna, min.edges = floor(n.loci*0.01), phi = 8, rm.COCs=FALSE) # for summary
clusters  <- extractClusters(ldna, min.edges = floor(n.loci*0.01), phi = 8, rm.COCs=TRUE)
dev.off()

# get summary. 
ld.sum <- summaryLDna(ldna, clusters0, mat)

# write clusters dataframe.
clstdf <- as.data.frame(stri_list2matrix(clusters$clusters))
colnames(clstdf)<-names(clusters$clusters)
clstdf <- clstdf[, names(sort(sapply(clusters$clusters, length), decreasing = T))] # reorder by n.loci
write.table(clstdf, file=paste(out, '.ldna_clusters', sep=''),
quote = F, row.names = F, sep="\t")

# output summaries.
write.table(file = paste(out, '.ldna_sum', sep=''), ld.sum, 
quote = F, row.names = F, sep="\t")

# output clusters.
for (c in seq(1,length(clusters$clusters))) {
	c.len  <- length(clusters$clusters[[c]])
	c.name <- names(clusters$clusters[c])
	ld     <- ld.sum[c.name, "Median.LD"]
    if (c.len >= floor(n.loci*0.01) && ld >= 0.3) {
        # out SNP positions
        c.pos.name <- paste(out, c.name, ld, c.len, sep='_')
        write.table(clusters$clusters[c], file=c.pos.name, quote = F, 
        row.names = F, col.names = F, sep="\t")
    }
}

# network pics
# pdf(paste(out,'.network.pdf', sep=''), height=7, width=5)
# par(mfcol=c(2,2))
# plotLDnetwork(ldna, mat, option=2, clusters=clusters[[1]], summary=ld.sum)
# dev.off()
