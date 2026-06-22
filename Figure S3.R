library(ggplot2)

setwd("C:/graduate/WH赤霉病/扩增子/WH Fungi")
rm(list=ls())
data<-read.csv("fun.relabu.Fusarium.csv",head = T)

p <- ggplot(data, aes(Sample,Fusarium)) +
  geom_col( aes(fill = Group),width = 0.9) +
  scale_fill_manual(values = c("#8b66b8","#c63596","#51c4c2","#2AB34A")) +
  labs(title = NULL, x = NULL, y = 'Relative abundance', fill = NULL) +
  ggtitle("Fusarium")+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 20,face="bold"),
        axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16,angle = -20),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of fusarium each sample.pdf", p,
       height=8,width =12,limitsize = FALSE )