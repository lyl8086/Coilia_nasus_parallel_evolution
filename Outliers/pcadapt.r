library(pcadapt)
library(qvalue)
filename <- read.pcadapt('input.bed', type = "bed")
res <- pcadapt(filename, K = 20, LD.clumping = list(size = 50, thr = 0.1))
plot(res, option = "screeplot")
dev.off()
k<-2
res <- pcadapt(filename, K = as.numeric(k), LD.clumping = list(size = 50, thr = 0.1))
qval <- qvalue(res$pvalues)$qvalues
outliers <- which(qval <= 0.1)

write.table(outliers, "pcadapt.outliers.id", col.names=F, row.names=F, quote=F)
