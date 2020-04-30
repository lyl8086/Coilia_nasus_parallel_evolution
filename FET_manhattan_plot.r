library('ggplot2')
library('magrittr')
library('dplyr')
library('gtools')
library(ggthemr)
ggthemr('fresh', type="inner")
P<-read.table("12W_FET.P",header=T, comment.char="*")

gwasResults<-P
names(gwasResults)[1:2] <- c("CHR", "BP")

don <<- gwasResults %>% 
  
  # Compute chromosome size
  group_by(CHR=factor(CHR, levels = c(paste(rep("LG", 24), seq(1,24,1), sep="")))) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


plt<-function(dat=don, P=P, xlab="chromosome", ylab="-log10(p value)",tt="") {

    dat<- dat %>% select(CHR,BPcum,P)
    names(dat)[3]<-"Pval"
    p<-ggplot(dat, aes(x=BPcum)) +
    
    # Show all points
    geom_point(aes(y=-log10(Pval), color=factor(CHR, levels = c(paste(rep("LG", 24), seq(1,24,1), sep="")))), alpha=0.8, size=0.3) +
    scale_color_manual(values = rep(c("grey", "black"), 24 )) +
    # custom X axis:
    scale_x_continuous(expand = c(0.02, 0.01), label = axisdf$CHR, breaks= axisdf$center) +
    scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis
  
    # Custom the theme:
    labs(x=xlab,y=ylab, title=tt)+
    geom_hline(aes(yintercept=4), colour="red", size=0.3)+
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
    if(xlab=="") {
        p<-p+ theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
    }
    return(p)
}
th<-plt(dat=don, P="TH", xlab="", tt="Taihu Lake")
ch<-plt(dat=don, P="CH", xlab="", tt="Chaohu Lake")
hz<-plt(dat=don, P="HZ", xlab="", tt="Hongze Lake")
lm<-plt(dat=don, P="LM", tt="Luoma Lake")

library(gridExtra)
tiff("all_FET.tiff", width=1800, height=900)
grid.arrange(th,ch,hz,lm,ncol=1)
dev.off()

pdf("all_FET.pdf", width=18, height=9,pointsize=0.1)
grid.arrange(th,ch,hz,lm,ncol=1)
dev.off()

