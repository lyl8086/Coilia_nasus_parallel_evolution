argv   <- commandArgs(TRUE)
infile <- argv[1]
out    <- argv[2]
vcf    <- argv[3]
keep   <- argv[4]
cores  <- argv[5]
min.e  <- 30
max.phi <- 10
min.ld <- 0.3
library(LDna)
library(stringi)
library(data.table)
if(missing(cores)|is.na(cores)) cores=1
cat("Using ", cores, " cpus\n")
cat(paste("Read file", infile, "..."))
x<-fread(infile, header = T, na.strings = c("-nan","NA"))
names(x)[5] <- 'R.2'
x$R.2[is.na(x$R.2)]<-0
cat("done.\n")

# init matrix
cat("Initiating matrix...")
mat.names <- unique(c(x$POS1,tail(x$POS2,1)))
n.loci <- length(mat.names)
mat<-matrix(nrow=n.loci, ncol=n.loci, dimnames=list(mat.names,mat.names))
# fill lower tri 
mat[lower.tri(mat)]<-x$R.2
cat("done.\n")

# read raw data.
cat("Convert raw matrix...")
ldna <- LDnaRaw(mat, digits=3, mc.cores=cores)
cat("\n")
pdf(paste(out,'LDna.Tree.pdf', sep=''), width=7, height=5)
par(mfcol=c(1,2))

# optimize phi start from 2.
pre_clsts <- 0
j <- 0
for (i in 2:max.phi) {
    cat("phi=", i)
    opt.clst  <- extractClusters(ldna, min.edges = min.e, phi = i, rm.COCs=TRUE)
    len <- length(opt.clst$clusters)
    cat(" size=", len, "\n")
	if (len == pre_clsts) j <<- j + 1
    if ( j == 2 | i == max.phi) {
        phi <<- i
        dev.off()
        break
    }
    pre_clsts <<- len
}
clusters <- extractClusters(ldna, min.edges = min.e, phi = phi, rm.COCs=TRUE)
clusters0 <- extractClusters(ldna, min.edges = min.e, phi = phi, rm.COCs=FALSE) # for summary

# get summary. 
ld.sum <- summaryLDna(ldna, clusters0, mat)

# if only one cluster.
n.clusters <- length(clusters$clusters)
if (n.clusters == 1) {
    names(clusters$clusters) <- rownames(ld.sum[which(ld.sum$Type=="SOC"),])
}
    
# write clusters dataframe.
clstdf <- as.data.frame(stri_list2matrix(clusters$clusters))
names(clstdf)<-names(clusters$clusters)
clstdf <- clstdf[, names(sort(sapply(clusters$clusters, length), decreasing = T)), drop=F] # reorder by n.loci
write.table(clstdf, file=paste(out, '.ldna_clusters', sep=''),
quote = F, row.names = F, sep="\t")

# output summaries.
write.table(file = paste(out, '.ldna_sum', sep=''), ld.sum, 
quote = F, row.names = F, sep="\t")

# output clusters.
for (c in 1:n.clusters) {
	c.len  <- length(clusters$clusters[[c]])
	c.name <- names(clusters$clusters[c])
	ld     <- ld.sum[c.name, "Median.LD"]
    if (c.len >= min.e & ld >= min.ld) {
        # out SNP positions
        c.pos <- list(chr=rep(x$CHR[1], c.len), pos=clusters$clusters[[c]])
        c.pos.name <- paste(out, c.name, ld, c.len, sep='_')
        write.table(c.pos, file=c.pos.name, quote = F, 
        row.names = F, col.names = F, sep="\t")
        # extract vcf
        system(paste("vcftools --gzvcf", vcf, 
        "--positions", c.pos.name,
        "--keep", keep,
        "--recode --out", c.pos.name, sep=" "))
    }
}
