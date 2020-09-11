library(adegenet)
library(factoextra)
library(dplyr)
library(ggsci)
dat<-read.genepop("input.gen", ncode=3)
#
levels(pop(dat))<-c("Yangtze River Estuary","Taihu Lake","Chaohu Lake","Hongze Lake","Luoma Lake")

dat.tab<-tab(dat, freq = TRUE, NA.method = "mean")
dat.pca<-dudi.pca(df = dat.tab, scale = FALSE, scannf = FALSE, nf = 2)
dat.mean<-dat.pca$li %>% 
  mutate(pop=pop(dat)) %>%
  group_by(pop) %>% 
  mutate(avg1=mean(Axis1), avg2=mean(Axis2),sd1=sd(Axis1),sd2=sd(Axis2))

dat.pcam<-dat.pca
dat.pcam$li$Axis1<-dat.mean$avg1
dat.pcam$li$Axis2<-dat.mean$avg2
p1<<-fviz_pca_ind(dat.pca, col.ind = pop(dat), 
                 addEllipses = TRUE, ellipse.level=0.667, 
                 repel = TRUE, 
                 legend.title = "Populations",
                 title="", label="")+scale_color_npg()+
				 theme(text=element_text(family="Helvetica", size=14),
				 legend.position="left")
p2<-fviz_pca_ind(dat.pcam, col.ind = pop(dat), 
                 xlim=c(min(dat.pca$li$Axis1), max(dat.pca$li$Axis1)),
                 ylim=c(min(dat.pca$li$Axis2), max(dat.pca$li$Axis2)),
                 repel = TRUE, 
                 legend.title = "Populations",
                 title="", label="")+scale_color_npg()
ggsave("ind_pca.pdf",p1,width=6, height=5)
ggsave("indm_pca.pdf",p2)
