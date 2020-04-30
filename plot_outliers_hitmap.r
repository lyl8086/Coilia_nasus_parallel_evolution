library(ggplot2)
library(reshape2)
library(dplyr)
library(ggthemr)
#ggthemr('fresh', type="inner")
library("RColorBrewer")
x<-read.table('FWA.frq', header=T, comment.char="*")
names(x)[1]<-'CHROM'

dat<-mutate(x, Average=(CH+TH+HZ+LM)/4)
dat<-melt(dat, id=c("CHROM", "POS", "REF", "ALT", "Major", "CJ"))
dat$variable<-factor(dat$variable, levels=c("TH","CH","HZ","LM","Average"), 
labels=c("Taihu Lake","Chaohu Lake","Hongze Lake","Luoma Lake","Average"))
#dat<-rbind(dat2, dat)

p<-ggplot(dat, aes(value,CJ))+ylim(-0.01,1.01)+xlim(-0.01,1.01)+
facet_wrap(~variable, ncol=1)+
#bins=60
geom_bin2d(binwidth=0.025)+
scale_fill_gradientn(colors=c(brewer.pal(9,"Blues")[2], rainbow(7)), breaks=c(1,50,100,250,300))+
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
ggsave('outliers_hitmap.pdf', p, width=4,height=10)

# histo

dat<-mutate(x, TH_Change=abs(CJ-TH), CH_Change=abs(CJ-CH), HZ_Change=abs(CJ-HZ), LM_Change=abs(CJ-LM), 
Average_Change=abs(CJ-(TH+CH+HZ+LM)/4)) %>%
select(-CJ,-TH,-CH,-HZ,-LM)

#dat<-melt(dat, id=c("CHROM", "POS", "REF", "ALT"))
#dat$variable<-factor(dat$variable, levels=c("TH_Change","CH_Change","HZ_Change","LM_Change"))
pop<-c("Taihu Lake","Chaohu Lake","Hongze Lake","Luoma Lake","Average")
names(pop)<-c('TH_Change', 'CH_Change', 'HZ_Change', 'LM_Change', 'Average_Change')

plt<- function(i,j) {
p<-ggplot(dat)+
geom_histogram(aes_string(x=i), breaks=seq(0,1,0.05), closed = "right")+
scale_x_continuous( label = seq(0,1,0.1), breaks=seq(0,1,0.1), limits=c(-0.01,1.01))+
labs(x="", y=j, title=pop[[i]])+
theme(
text=element_text(family="Helvetica"),
plot.title = element_text(hjust = 0.5,size=rel(1)),
panel.background = element_blank(),
panel.border = element_blank(),
axis.line = element_line(colour = "black"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank()
)
return(p)
}

p1<-plt('TH_Change', "count")
p2<-plt('CH_Change', "count")
p3<-plt('HZ_Change', "count")
p4<-plt('LM_Change', "count")
p<-plt('Average_Change', "count")
library(gridExtra)
pdf('allele_change.pdf', height=6,width=3)
#grid.arrange(grobs=lapply(names(dat)[6:ncol(dat)],myplot), ncol=2)
grid.arrange(p1,p2,p3,p4,ncol=1)
dev.off()

#ggsave('outliers_hist.pdf', p, width=4,height=10)
