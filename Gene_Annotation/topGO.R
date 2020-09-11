library(topGO)

# prepare data.
geneID2GO <- readMappings("genes_to_go")
geneNames <- names(geneID2GO)
myInterestingGenes <- read.table("gene_of_interst", header=F, stringsAsFactors=F)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes$V1))
names(geneList) <- geneNames

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
annot = annFUN.gene2GO, gene2GO = geneID2GO)
allGO <- usedGO(object = GOdata) 

# do enrichment.
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher_w <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
# resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
# resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(GOdata, topNodes = length(allGO), numChar=1000, 
	classicFisher = resultFisher,
	weightFisher = resultFisher_w,
	#classicKS = resultKS, 
	#elimKS = resultKS.elim,
	orderBy = "weightFisher", ranksOf = "classicFisher")
fdr <- round(p.adjust(allRes$weightFisher, method="BH"), digits = 5)
allRes <- cbind(allRes, fdr)
# write the result.
write.table(allRes, file="GO_enrichment.tsv", quote=F, row.names=F, col.names=T, sep="\t")

