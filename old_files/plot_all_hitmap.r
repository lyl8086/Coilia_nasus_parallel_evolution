library(ggplot2)
library(reshape2)
library(dplyr)
library("RColorBrewer")
library(ggthemr)
x<-read.table('all.maf', header=T, comment.char="*")
names(x)[1]<-'CHROM'
dat<-mutate(x, Average=(CH+TH+HZ+LM)/4)
dat<-melt(dat, id=c("CHROM", "POS", "REF", "ALT", "Major", "CJ"))
dat$variable<-factor(dat$variable, levels=c("TH","CH","HZ","LM","Average"), 
labels=c("Taihu Lake","Chaohu Lake","Hongze Lake","Luoma Lake","Average"))

p<-ggplot(dat, aes(value,CJ))+ylim(0.49,1.01)+xlim(-0.01,1.01)+
facet_wrap(~variable, ncol=1)+
geom_bin2d(binwidth=0.025)+
scale_fill_gradientn(colors=c(brewer.pal(9,"Blues")[2],rainbow(7)), 
breaks=c(1,500,1000,2000,2500,3000,4000,5000), 
limits=c(1,5500), na.value=rainbow(8)[8])+
labs(x="", y="Yangtze River Estuary")+
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
  
#ggthemr('fresh')
ggsave('all_hitmap.pdf', p, width=4,height=10)
