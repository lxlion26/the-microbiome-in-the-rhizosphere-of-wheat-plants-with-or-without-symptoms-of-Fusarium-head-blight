#################Bacterial beta dispersion###################
library(vegan)
library(plotrix)
library(ggbreak)
library(ade4)
library(adespatial)
library(ggtern)

setwd("C:/graduate/WH赤霉病/扩增子/WH Bacteria")
rm(list=ls())

otu <- read.csv("da.rar.csv",head=T,row.names=1)
env<-read.csv("env.csv",head=T,row.names=1)
groups <- factor(env$State)

dis <- vegdist(otu, method = 'bray')
mod <- betadisper(d = dis, group = groups, type = 'centroid')
shapiro.test(mod$distances)#×
anova(mod)
set.seed(123)
permutest(mod, pairwise = TRUE, permutations = 999)

#####PCoA#######
plot(mod,ellipse = TRUE, hull = FALSE, conf = 0.95,main = NULL, sub = NULL,col = c("#8b66b8","#51c4c2"))
grid(col = "lightgray")

boxplot(mod,col = c("#8b66b8","#51c4c2"), ylab = 'Beta Dispersion',xlab=NULL)

mod<-as.data.frame(mod$distances)
mod<-cbind(env$Sample,mod[,1],env$Compartment,env$State)
mod<-as.data.frame(mod)
colnames(mod)=c("Sample","distances","Compartment","State")
str(mod)
mod$distances<-as.numeric(mod$distances)
str(mod)
write.csv(mod, 'bacteria.betadisper.csv', quote = FALSE)

otu1<-otu[env$Compartment=="Rhizosphere",]
env1<-env[env$Compartment=="Rhizosphere",]
groups <- factor(env1$State)

dis <- vegdist(otu1, method = 'bray')

mod1 <- betadisper(d = dis, group = groups, type = 'centroid')

shapiro.test(mod1$distances)#×
anova(mod1)

set.seed(123)
permutest(mod1, pairwise = TRUE, permutations = 999)

otu2<-otu[env$Compartment=="Root",]
env2<-env[env$Compartment=="Root",]
groups <- factor(env2$State)

dis <- vegdist(otu2, method = 'bray')

mod2 <- betadisper(d = dis, group = groups, type = 'centroid')

shapiro.test(mod2$distances)#对
anova(mod2)

set.seed(123)
permutest(mod2, pairwise = TRUE, permutations = 999)

p<-ggplot(mod, aes(x = State, y = distances,color=State)) +
  geom_boxplot(size=0.8, width = 0.3,alpha=0) +
  geom_jitter(position =position_jitter(0.4),alpha=0.2,size=1.5)+
  scale_color_manual(values=c("#8b66b8","#51c4c2"))+
  labs(title = "", x = NULL, y = 'Distance to centroid', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  scale_y_continuous(expand = c(0, 0), limit = c(0.3, 0.85))+
  ggtitle("Bacterial community",expression('Rhizosphere:'~italic(F)~'= 3.240'^ns~~~~~'Root:'~italic(F) ~'= 1.775'^ns))+
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
        plot.title=element_text(face = "bold", colour = "black", size = 20),
        plot.subtitle=element_text(face = "bold", colour = "black", size = 20))+
  facet_grid( .~env$Compartment)
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S beta dispersion.pdf", p,
       height=8,width =9,limitsize = FALSE )

#########################Fungal beta dispersion########################
setwd("C:/graduate/WH赤霉病/扩增子/WH Fungi")
rm(list=ls())

otu <- read.csv("da.rar.csv",head=T,row.names=1)
env<-read.csv("env.csv",head=T,row.names=1)
groups <- factor(env$State)

dis <- vegdist(otu, method = 'bray')
mod <- betadisper(d = dis, group = groups, type = 'centroid')
shapiro.test(mod$distances)#×
anova(mod)

set.seed(123)
permutest(mod, pairwise = TRUE, permutations = 999)

#####PCoA#######
plot(mod,ellipse = TRUE, hull = FALSE, conf = 0.95,main = NULL, sub = NULL,col = c("#8b66b8","#51c4c2"))
grid(col = "lightgray")

boxplot(mod,col = c("#8b66b8","#51c4c2"), ylab = 'Beta Dispersion',xlab=NULL)

mod<-as.data.frame(mod$distances)
mod<-cbind(env$Sample,mod[,1],env$Compartment,env$State)
mod<-as.data.frame(mod)
colnames(mod)=c("Sample","distances","Compartment","State")
str(mod)
mod$distances<-as.numeric(mod$distances)
str(mod)
write.csv(mod, 'fungi.betadisper.csv', quote = FALSE)

otu1<-otu[env$Compartment=="Rhizosphere",]
env1<-env[env$Compartment=="Rhizosphere",]
groups <- factor(env1$State)

dis <- vegdist(otu1, method = 'bray')

mod1 <- betadisper(d = dis, group = groups, type = 'centroid')

shapiro.test(mod1$distances)#×
anova(mod1)


set.seed(123)
permutest(mod1, pairwise = TRUE, permutations = 999)

otu2<-otu[env$Compartment=="Root",]
env2<-env[env$Compartment=="Root",]
groups <- factor(env2$State)
dis <- vegdist(otu2, method = 'bray')

mod2 <- betadisper(d = dis, group = groups, type = 'centroid')

shapiro.test(mod2$distances)#×
anova(mod2)

set.seed(123)
permutest(mod2, pairwise = TRUE, permutations = 999)

p<-ggplot(mod, aes(x = State, y = distances,color=State)) +
  geom_boxplot(size=0.8, width = 0.3,alpha=0) +
  geom_jitter(position =position_jitter(0.4),alpha=0.2,size=1.5)+
  scale_color_manual(values=c("#8b66b8","#51c4c2"))+
  labs(title = "", x = NULL, y = 'Distance to centroid', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  scale_y_continuous(expand = c(0, 0), limit = c(0.2, 0.8))+
  ggtitle("Fungal community",expression('Rhizosphere:'~italic(F)~'= 0.023'^ns~~~~~'Root:'~italic(F) ~'= 0.193'^ns))+
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
        plot.title=element_text(face = "bold", colour = "black", size = 20),
        plot.subtitle=element_text(face = "bold", colour = "black", size = 20))+
  facet_grid( .~env$Compartment)
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS beta dispersion.pdf", p,
       height=8,width =9,limitsize = FALSE )