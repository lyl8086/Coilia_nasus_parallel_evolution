#!/usr/bin/R

library(lostruct)
library(colorspace)
library(jsonlite)
library(RColorBrewer)
library(gtools)
pop_pairs = c("CJ+TH", "CJ+CH", "CJ+HZ", "CJ+LM")
pop_main  = c(
"Taihu Lake vs. Yangtze River Estuary",
"Chaohu Lake vs. Yangtze River Estuary",
"Hongze Lake vs. Yangtze River Estuary",
"Luoma Lake vs. Yangtze River Estuary")
pdf("local_PCA.MDS1.pdf", width=17, height=9)
par(mfrow=c(4,1))
#' Set up to plot all chromosomes together
#' will plot the vector of values 'y' along all the chromosomes.
chrom.plot <- function (y,ylab='',main='',chrom.labels=TRUE,...) {
    plot(0, type='n', xlim=range(chrom.offsets/1e6), ylim=range(y,finite=TRUE), 
         xlab='', xaxt='n', yaxt='n', ylab=ylab, main=main)
    if (length(chroms)>1) for (k in 1:floor(length(chroms)/2)) {
        rect( xleft=chrom.dividers[2*k-1]/1e6, xright=chrom.dividers[2*k]/1e6, 
             ybottom=par("usr")[3], ytop=par("usr")[4], 
             border=NA, col=adjustcolor("grey",0.25) )
    }
    abline( v=chrom.dividers/1e6, lty=3, col=adjustcolor("grey",0.5), lwd=2 )
    if (chrom.labels) axis( 1, at=chrom.mids/1e6, labels=chroms, las=0, tick=FALSE )
    points( regions$pos/1e6, y, ...)
}

for (i in 1:length(pop_pairs)) {
cur_path <- paste0(pop_pairs[i],".out/res")
# data files, precomputed
mds.file <- paste0(cur_path, "/", "mds_coords.csv")
regions.files <- mixedsort(list.files(paste0(cur_path, "/"),".*.regions.csv"))

# read in mds
mds.coords <- read.csv(mds.file, header=TRUE)
mds.cols <- (1:ncol(mds.coords))[-(1:2)]

# sort chrom.
mds.coords <- mds.coords[mixedorder(as.character(mds.coords$chrom)),]

# position information
regions <- do.call( rbind, lapply(paste0(cur_path, "/", regions.files), read.csv, header=TRUE ) )
# figure out where to plot things at
chroms <- unique(regions$chrom)
chrom.starts <- tapply( regions$start, regions$chrom, min, na.rm=TRUE )
chrom.ends <- tapply( regions$end, regions$chrom, max, na.rm=TRUE )
chrom.spacing <- floor(.05*mean(chrom.ends))
chrom.offsets <- c(0,cumsum(chrom.spacing+chrom.ends))
names(chrom.offsets) <- c(chroms,"end")
chrom.dividers <- c(0,chrom.offsets[-1])-chrom.spacing/2
chrom.mids <- chrom.dividers[-1] - diff(chrom.dividers)/2

# this is where to plot windows at when plotting with all chromosomes
regions$pos <- chrom.offsets[match(regions$chrom,chroms)]+(regions$start+regions$end)/2

chrom.cols <- rainbow_hcl(length(chroms), c=90, end=.9*360)[as.numeric(regions$chrom)]

#Here are the extreme windows in the MDS plot:
mds.corners <- corners( mds.coords[,mds.cols[1:2]], prop=.05, k=1)
# set up colors and pchs for corners
corner.cols <- brewer.pal(3,"Dark2")
corner.pch <- c(15,17,19)
ccols <- rep("black",nrow(mds.coords))
cpch <- rep(20,nrow(mds.coords))
for (k in 1:ncol(mds.corners)) {
    ccols[ mds.corners[,k] ] <- corner.cols[2]
    cpch[ mds.corners[,k] ] <- corner.pch[k]
}

# plot the first mds
k <- mds.cols[1]
chrom.plot(		
		mds.coords[,k], 
		main=pop_main[i],pch=20, 
        xlab="Position (Mb)",
        chrom.labels=TRUE,
        ylab=colnames(mds.coords)[k],
        col=adjustcolor(ccols,0.75) )
}
