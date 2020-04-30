library('ggplot2')
library('magrittr')
library('dplyr')
library('gtools')
library('reshape2')
library(gridExtra)
library(ggthemr)
#ggthemr('fresh')
x<-read.table('FWA.frq', header=T, comment.char="+")

names(x)<-c('CHROM','POS','REF','ALT','FWA',"CJ", "CH","HZ","LM","TH")

plt<-function(i="",j="") {

    p<-ggplot(x, aes_string(i))+
    geom_histogram(breaks=seq(0,1,0.05))+
    scale_x_continuous( label = seq(0,1,0.1), breaks=seq(0,1,0.1), limits=c(0,1))+
    labs(title=j, x="")+
    theme(
text=element_text(family="Helvetica"),
plot.title = element_text(hjust = 0.5,size=rel(1)),
panel.background = element_blank(),
panel.border = element_blank(),
axis.line = element_line(colour = "black"),
panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(),
strip.background = element_blank()
)  
    return(p)
}
p<-plt('CJ', "Yangtze River Estuary")
p1<-plt('TH', "Taihu Lake")
p2<-plt('CH', "Chaohu Lake")
p3<-plt('HZ', "Hongze Lake")
p4<-plt('LM', "Luoma Lake")
pdf('fwa_1269.pdf', height=8,width=3)
grid.arrange(p,p1,p2,p3,p4,ncol=1)
dev.off()
