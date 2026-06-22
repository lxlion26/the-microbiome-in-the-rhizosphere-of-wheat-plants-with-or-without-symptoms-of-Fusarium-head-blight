####################Bacterial network properties#########################
library(vegan)
library(ggplot2)
library(colorRamps)
library(ape)
library(splitstackshape)
library(reshape2)

setwd("C:/graduate/WH赤霉病/扩增子/WH Bacteria")
rm(list=ls())
BBNetwork <- read.csv("BB Network（正负）.csv", head = T, row.names = 1)

p<-ggplot(BBNetwork, aes(x = Compartment, y = connectance, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=connectance),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Bacterial Network", x = NULL, y = 'Connectance', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 0.020))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S network connectance.pdf", p,
       height=8,width =9,limitsize = FALSE )

p<-ggplot(BBNetwork, aes(x = Compartment, y = average.degree, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=average.degree),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Bacterial Network", x = NULL, y = 'Average degree', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 4))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S network average.degree.pdf", p,
       height=8,width =9,limitsize = FALSE )

p<-ggplot(BBNetwork, aes(x = Compartment, y = average.path.length, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=average.path.length),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Bacterial Network", x = NULL, y = 'Average path length', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 6))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S network average.path.length.pdf", p,
       height=8,width =9,limitsize = FALSE )

p<-ggplot(BBNetwork, aes(x = Compartment, y = diameter, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=diameter),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Bacterial Network", x = NULL, y = 'Diameter', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 15))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S network diameter.pdf", p,
       height=8,width =9,limitsize = FALSE )

p<-ggplot(BBNetwork, aes(x = Compartment, y = clustering.coefficient, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=clustering.coefficient),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Bacterial Network", x = NULL, y = 'Clustering coefficient', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 0.5))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S network clustering.coefficient.pdf", p,
       height=8,width =9,limitsize = FALSE )

p<-ggplot(BBNetwork, aes(x = Compartment, y = no.clusters, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=no.clusters),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Bacterial Network", x = NULL, y = 'Number of clusters', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 80))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S network no.clusters.pdf", p,
       height=8,width =9,limitsize = FALSE )

p<-ggplot(BBNetwork, aes(x = Compartment, y = centralization.betweenness, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=centralization.betweenness),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Bacterial Network", x = NULL, y = 'Centralization betweenness', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 0.12))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S network centralization.betweenness.pdf", p,
       height=8,width =9,limitsize = FALSE )

p<-ggplot(BBNetwork, aes(x = Compartment, y = centralization.degree, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=centralization.degree),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Bacterial Network", x = NULL, y = 'Centralization degree', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 0.07))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S network centralization.degree.pdf", p,
       height=8,width =9,limitsize = FALSE )

p<-ggplot(BBNetwork, aes(x = Compartment, y = Modularity, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=Modularity),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Bacterial Network", x = NULL, y = 'Modularity', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 1.03)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 1))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S network Modularity.pdf", p,
       height=8,width =9,limitsize = FALSE )

p<-ggplot(BBNetwork, aes(x = Compartment, y = No.modules, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=No.modules),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Bacterial Network", x = NULL, y = 'Number of modules', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 80))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S network No.modules.pdf", p,
       height=8,width =9,limitsize = FALSE )

########################Fungal network properties#####################
setwd("C:/graduate/WH赤霉病/扩增子/WH Fungi")

rm(list=ls())
FFNetwork <- read.csv("FF Network（正负）.csv", head = T, row.names = 1)

p<-ggplot(FFNetwork, aes(x = Compartment, y = connectance, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=connectance),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Fungal Network", x = NULL, y = 'Connectance', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 0.06))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS network connectance.pdf", p,
       height=8,width =9,limitsize = FALSE )

p<-ggplot(FFNetwork, aes(x = Compartment, y = average.degree, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=average.degree),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Fungal Network", x = NULL, y = 'Average degree', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 5))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS network average.degree.pdf", p,
       height=8,width =9,limitsize = FALSE )

p<-ggplot(FFNetwork, aes(x = Compartment, y = average.path.length, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=average.path.length),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Fungal Network", x = NULL, y = 'Average path length', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 3))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS network average.path.length.pdf", p,
       height=8,width =9,limitsize = FALSE )

p<-ggplot(FFNetwork, aes(x = Compartment, y = diameter, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=diameter),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Fungal Network", x = NULL, y = 'Diameter', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 10))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS network diameter.pdf", p,
       height=8,width =9,limitsize = FALSE)

p<-ggplot(FFNetwork, aes(x = Compartment, y = clustering.coefficient, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=clustering.coefficient),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Fungal Network", x = NULL, y = 'Clustering coefficient', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 1))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS network clustering.coefficient.pdf", p,
       height=8,width =9,limitsize = FALSE)

p<-ggplot(FFNetwork, aes(x = Compartment, y = no.clusters, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=no.clusters),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Fungal Network", x = NULL, y = 'Number of clusters', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 40))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS network no.clusters.pdf", p,
       height=8,width =9,limitsize = FALSE)

p<-ggplot(FFNetwork, aes(x = Compartment, y = centralization.betweenness, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=centralization.betweenness),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Fungal Network", x = NULL, y = 'Centralization betweenness', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 0.02))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS network centralization.betweenness.pdf", p,
       height=8,width =9,limitsize = FALSE)

p<-ggplot(FFNetwork, aes(x = Compartment, y = centralization.degree, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=centralization.degree),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Fungal Network", x = NULL, y = 'Centralization degree', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 0.2))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS network centralization.degree.pdf", p,
       height=8,width =9,limitsize = FALSE)

p<-ggplot(FFNetwork, aes(x = Compartment, y = Modularity, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=Modularity),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Fungal Network", x = NULL, y = 'Modularity', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 1.03)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 1))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS network Modularity.pdf", p,
       height=8,width =9,limitsize = FALSE)

p<-ggplot(FFNetwork, aes(x = Compartment, y = No.modules, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=No.modules),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Fungal Network", x = NULL, y = 'Number of modules', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 40))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS network No.modules.pdf", p,
       height=8,width =9,limitsize = FALSE)
