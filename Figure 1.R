library(picante)
library(ggplot2)
library(ggbreak)

####################Bacterial Shannon#######################################
rm(list=ls())
setwd("C:/graduate/WH赤霉病/扩增子/WH Bacteria")

otu <- read.csv("da.rar.csv",head=T,row.names=1)
env<-read.csv("env.csv",head=T,row.names=1)
metadata.otu<-cbind(env,otu)

alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}
alpha_all <- alpha(otu, base = 2)

write.csv(alpha_all, 'bacteria.alpha.csv', quote = FALSE)
alpha <- read.csv('bacteria.alpha.csv',head=T,row.names=1)
alpha<-cbind(env,alpha)

Rhizosphere<-alpha[alpha$Compartment=='Rhizosphere',]
shapiro.test(Rhizosphere$Shannon) #×
shapiro.test(log(Rhizosphere$Shannon)) #×
shapiro.test(sqrt(Rhizosphere$Shannon))# ×

kruskal.test(Rhizosphere$Shannon ~ Rhizosphere$State)#

Root<-alpha[alpha$Compartment=='Root',]
shapiro.test(Root$Shannon) #×
shapiro.test(log(Root$Shannon)) #×
shapiro.test(sqrt(Root$Shannon))# ×

kruskal.test(Root$Shannon ~ Root$State)#

p<-ggplot(alpha, aes(x = State, y = Shannon,color=State)) +
  geom_boxplot(size=0.8, width = 0.3,alpha=0) +
  geom_jitter(position =position_jitter(0.3),alpha=0.2,size=1.5)+
  scale_color_manual(values=c("#8b66b8","#51c4c2"))+
  labs(title = "", x = NULL, y = 'Bacterial Shannon diversity index', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  scale_y_continuous(expand = c(0, 0), limit = c(1, 9.5),breaks=seq(1,9.5,1))+
  ggtitle(expression(~~~~'Rhizosphere:'~italic(χ) ^2~'= 4.547*'~~~~~'Root:'~italic(χ) ^2~'= 0.026'^ns))+
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
        plot.title=element_text(face = "bold", colour = "black", size = 20))+
  facet_grid( .~env$Compartment)+
  scale_y_break(breaks = c(2, 3), scales = 5, ticklabels=seq(3,9.5,1))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S shannon.pdf", p,
       height=8,width =9,limitsize = FALSE )

########################Fungal Shannon################################################
rm(list=ls())
setwd("C:/graduate/WH赤霉病/扩增子/WH Fungi")

otu <- read.csv("da.rar.csv",head=T,row.names=1)
env<-read.csv("env.csv",head=T,row.names=1)
metadata.otu<-cbind(env,otu)

alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}
alpha_all <- alpha(otu, base = 2)

write.csv(alpha_all, 'fungi.alpha.csv', quote = FALSE)
alpha <- read.csv('fungi.alpha.csv',head=T,row.names=1)
alpha<-cbind(env,alpha)

Rhizosphere<-alpha[alpha$Compartment=='Rhizosphere',]
shapiro.test(Rhizosphere$Shannon) #×
shapiro.test(log(Rhizosphere$Shannon)) #×
shapiro.test(sqrt(Rhizosphere$Shannon))# ×

kruskal.test(Rhizosphere$Shannon ~ Rhizosphere$State)#

Root<-alpha[alpha$Compartment=='Root',]
shapiro.test(Root$Shannon) #×
shapiro.test(log(Root$Shannon)) #×
shapiro.test(sqrt(Root$Shannon))# ×

kruskal.test(Root$Shannon ~ Root$State)#

p<-ggplot(alpha, aes(x = State, y = Shannon,color=State)) +
  geom_boxplot(size=0.8, width = 0.3,alpha=0) +
  geom_jitter(position =position_jitter(0.3),alpha=0.2,size=1.5)+
  scale_color_manual(values=c("#8b66b8","#51c4c2"))+
  labs(title = "", x = NULL, y = 'Fungal Shannon diversity index', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  scale_y_continuous(expand = c(0, 0), limit = c(2, 8))+
  ggtitle(expression(~~'Rhizosphere:'~italic(χ) ^2~'= 3.686'^ns~~~~~'Root:'~italic(χ) ^2~'= 1.060'^ns))+
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
        plot.title=element_text(face = "bold", colour = "black", size = 20))+
  facet_grid( .~env$Compartment)
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS shannon.pdf", p,
       height=8,width =9,limitsize = FALSE )

######################Bacterial CCA###############################
library(vegan)
library(ggplot2)
library(ape)
library(splitstackshape)
library(reshape2)
library(dplyr)

setwd("C:/graduate/WH赤霉病/扩增子/WH Bacteria")
rm(list = ls())

da.rar<-read.csv("da.rar.csv",head = T, row.names = 1)
env<-read.csv("env.csv",head = T, row.names = 1)
adonis2(da.rar~env$State * env$Compartment) 

bac.cca<-cca(da.rar~env$State * env$Compartment)
anova(bac.cca, by="term", permutations=499)
cca.tmp<-summary(bac.cca)

env.cca<-data.frame(cca.tmp$sites, env)

cca1<-round(bac.cca$CCA$eig[1],3)*100
cca2<-round(bac.cca$CCA$eig[2],3)*100

env.cca$Group<-paste(env.cca$State,env.cca$Compartment)

p<-ggplot(data=env.cca, aes(x= CCA1, y= CCA2 ))+ 
  geom_point(aes(colour= Group, shape = Compartment),size=4, alpha = 0.8) +
  labs(x = sprintf("CCA1 (%.1f%%)", cca1), y = sprintf("CCA2 (%.1f%%)", cca2))+
  scale_fill_manual(values= c("#8b66b8","#c63596","#51c4c2","#2AB34A"))+
  scale_colour_manual(values= c("#8b66b8","#c63596","#51c4c2","#2AB34A"))+
  stat_ellipse(aes(fill = Group), geom = 'polygon', level = 0.95, alpha = 0.1,show.legend = FALSE)+
  theme_bw() +
  ggtitle("Bacterial community",expression(italic(R) ["State"]^2~'= 0.009*'~~~~~italic(R)["Compartment"] ^2~'= 0.192***'))+
  theme(legend.title = element_text(colour="black", size=20, face="bold"),
        legend.text = element_text(colour="black", size=16),
        strip.text = element_text(size = 20,face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),legend.box="vertical",
        legend.position = "right",
        plot.title=element_text(face = "bold", colour = "black", size = 20),
        plot.subtitle=element_text(face = "bold", colour = "black", size = 20))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S CCA.pdf", p,
       height=9,width =12,limitsize = FALSE )

######################Fungal CCA###############################
setwd("C:/graduate/WH赤霉病/扩增子/WH Fungi")
rm(list = ls())

da.rar<-read.csv("da.rar.csv",head = T, row.names = 1)
env<-read.csv("env.csv",head = T, row.names = 1)
adonis2(da.rar~env$State * env$Compartment)  

fun.cca<-cca(da.rar~env$State * env$Compartment)
anova(fun.cca, by="term", permutations=499)
cca.tmp<-summary(fun.cca)

env.cca<-data.frame(cca.tmp$sites, env)

cca1<-round(fun.cca$CCA$eig[1],3)*100
cca2<-round(fun.cca$CCA$eig[2],3)*100

env.cca$Group<-paste(env.cca$State,env.cca$Compartment)

p<-ggplot(data=env.cca, aes(x= CCA1, y= CCA2 ))+ 
  geom_point(aes(colour= Group, shape = Compartment),size=4, alpha = 0.8) +
  labs(x = sprintf("CCA1 (%.1f%%)", cca1), y = sprintf("CCA2 (%.1f%%)", cca2))+
  scale_fill_manual(values= c("#8b66b8","#c63596","#51c4c2","#2AB34A"))+
  scale_colour_manual(values= c("#8b66b8","#c63596","#51c4c2","#2AB34A"))+
  stat_ellipse(aes(fill = Group), geom = 'polygon', level = 0.95, alpha = 0.1,show.legend = FALSE)+
  theme_bw() +
  ggtitle("Fungal community",expression(italic(R) ["State"]^2~'= 0.016***'~~~~~italic(R)["Compartment"] ^2~'= 0.076***'))+
  theme(legend.title = element_text(colour="black", size=20, face="bold"),
        legend.text = element_text(colour="black", size=16),
        strip.text = element_text(size = 20,face="bold"),
        panel.spacing = unit(0, "lines"),
        axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),legend.box="vertical",
        legend.position = "right",
        plot.title=element_text(face = "bold", colour = "black", size = 20),
        plot.subtitle=element_text(face = "bold", colour = "black", size = 20))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS CCA.pdf", p,
       height=9,width =12,limitsize = FALSE )

#####################PERMANOVA R2##################
library(vegan)
library(ggplot2)
library(colorRamps)
library(ape)
library(splitstackshape)
library(reshape2)

setwd("C:/graduate/WH赤霉病/扩增子")
rm(list=ls())
Permanova <- read.csv("Permanova R方.csv")
col11<-c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#97D1A0", "#43A743", "grey","#ff00ff","#00ff00", "deepskyblue", "gold", "red", "navy", "darkgreen","maroon3", "black", "bisque", "grey")
p<-ggplot(Permanova, aes(x = Community, y = PERMANOVA.R2, fill=Source.of.variation)) +
  geom_bar(stat='identity', position = "stack", width=0.5)+
  scale_fill_manual(values= col11)+
  theme_bw()+
  labs(x="",y="")+
  ylab(bquote(PERMANOVA~R^2))+
  guides(fill=guide_legend(title= "Source of variation"))+
  #scale_y_continuous(limit = c(0, 1))+
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 20,face="bold"),
        legend.title = element_text(colour="black", size=20, face="bold"),
        legend.text = element_text(colour="black", size=16),
        axis.text.y=element_text(colour="black",size=16,face = "bold"),
        axis.text.x=element_text(colour="black",size=16,face = "bold",angle = -20),
        axis.title=element_text(colour="black",size=20,face="bold"))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/PERMANOVA R2.pdf", p,
       height=8,width =9,limitsize = FALSE )

#################Relative abundance of Fusarium###################
library(vegan)
library(ggplot2)
library(colorRamps)
library(ape)
library(splitstackshape)
library(reshape2)

setwd("C:/graduate/WH赤霉病/扩增子/WH Fungi")
rm(list = ls())
ID <- read.csv("ID.onlyFun.csv", head = T, row.names = 1)
ID$OTU<-rownames(ID)
data0 <- read.csv("da.rar.csv", head = T, row.names = 1)
env <- read.csv("env.csv", head = T)


ID$OTU.ID<-paste(rownames(ID),ID$phylum, ID$genus, sep="_")
sum(env$Sample != rownames(data0))  # 0 说明env和Otu table 中sample编号顺序是一致的
sum(rownames(ID) != colnames(data0)) # 0 说明ID和Otu table 中OTU编号顺序是一致的

Flev<-ID[,"genus"] 
fung.lev<-data.frame(aggregate(t(data0),by=list(Flev) , sum))
rownames(fung.lev)<-fung.lev[,1]; 
fung.lev<-fung.lev[,-1]
data1<-data.frame(t(fung.lev))

total<-apply(data1, 1, sum); 
fung.relabu<-data.frame(lapply(data1, function(x) {  x / total  }) )

fung.RZ<-fung.relabu
env.RZ<-env
fung.RZ<-cbind(fung.RZ,env.RZ[,c(3,4)])
fung.R<-fung.RZ[fung.RZ$Compartment=="Rhizosphere",]
fung.Z<-fung.RZ[fung.RZ$Compartment=="Root",]


shapiro.test(fung.R$Fusarium)#×
shapiro.test(log(fung.R$Fusarium)) #×
shapiro.test(sqrt(fung.R$Fusarium))# ×
shapiro.test(fung.Z$Fusarium)#×
shapiro.test(log(fung.Z$Fusarium)) #×
shapiro.test(sqrt(fung.Z$Fusarium))#×

kruskal.test(fung.R$Fusarium ~ fung.R$State)#
kruskal.test(fung.Z$Fusarium ~ fung.Z$State)#

p<-ggplot(fung.RZ, aes(x = State, y = Fusarium, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Fusarium",expression('Rhizosphere:'~italic(χ) ^2~'=  0.093'^ns~~~~~'Root:'~italic(χ) ^2~'= 0.014'^ns))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Fusarium.pdf", p,
       height=8,width =9,limitsize = FALSE )

###############Bacterial AVD###############################
library(ggplot2)

setwd("C:/graduate/WH赤霉病/扩增子/WH Bacteria")
rm(list=ls())

otu <- read.csv("da.rar.without0.csv",head=T,row.names=1)
env<-read.csv("env.csv",head=T,row.names=1)

ai <- abs(otu-apply(otu, 1, mean))/apply(otu, 1, sd)
avd <- colSums(ai)/(1*nrow(otu))
env$AVD <- avd

env.R<-env[env$Compartment=="Rhizosphere",]

shapiro.test(env.R$AVD) #对
shapiro.test(log(env.R$AVD)) #×
shapiro.test(sqrt(env.R$AVD))# ×

bartlett.test(env.R$AVD ~ interaction(State),data=env.R)#对

anova_env.R<- aov(env.R$AVD ~ State,data = env.R)
summary(anova_env.R)#6.707*

env.Z<-env[env$Compartment=="Root",]

shapiro.test(env.Z$AVD) #×
shapiro.test(log(env.Z$AVD)) #×
shapiro.test(sqrt(env.Z$AVD))# ×

kruskal.test(env.Z$AVD ~ env.Z$State)#0.256ns

p<-ggplot(env, aes(x = State, y = AVD,color=State)) +
  geom_boxplot(size=0.8, width = 0.3,alpha=0) +
  geom_jitter(position =position_jitter(0.2),alpha=0.2,size=1.5)+
  scale_color_manual(values=c("#8b66b8","#51c4c2"))+
  labs(title = "", x = NULL, y = 'AVD', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  ggtitle("Bacterial community",expression('Rhizosphere:'~italic(F)~'= 6.707*'~~~~~'Root:'~italic(χ) ^2~'= 0.256'^ns))+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 20,face="bold"),
        axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16,angle=-20),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20),
        plot.subtitle=element_text(face = "bold", colour = "black", size = 20))+
  facet_grid( .~env$Compartment)
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S AVD.pdf", p,
       height=8,width =9,limitsize = FALSE )

################Fungal AVD#########################
library(ggplot2)

setwd("C:/graduate/WH赤霉病/扩增子/WH Fungi")
rm(list=ls())

otu <- read.csv("da.rar.without0.csv",head=T,row.names=1)
env<-read.csv("env.csv",head=T,row.names=1)

ai <- abs(otu-apply(otu, 1, mean))/apply(otu, 1, sd)
avd <- colSums(ai)/(1*nrow(otu))
env$AVD <- avd

env.R<-env[env$Compartment=="Rhizosphere",]

shapiro.test(env.R$AVD) #×
shapiro.test(log(env.R$AVD)) #×
shapiro.test(sqrt(env.R$AVD))# ×

kruskal.test(env.R$AVD ~ env.R$State)#2.148ns

env.Z<-env[env$Compartment=="Root",]

shapiro.test(env.Z$AVD) #对
shapiro.test(log(env.Z$AVD)) #对
shapiro.test(sqrt(env.Z$AVD))# 对

bartlett.test(env.Z$AVD ~ interaction(State),data=env.Z)#对

anova_env.Z<- aov(env.Z$AVD ~ State,data = env.Z)
summary(anova_env.Z)#0.050ns

###作图
p<-ggplot(env, aes(x = State, y = AVD,color=State)) +
  geom_boxplot(size=0.8, width = 0.3,alpha=0) +
  geom_jitter(position =position_jitter(0.2),alpha=0.2,size=1.5)+
  scale_color_manual(values=c("#8b66b8","#51c4c2"))+
  labs(title = "", x = NULL, y = 'AVD', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  ggtitle("Fungal community",expression('Rhizosphere:'~italic(χ)^2~'= 2.148'^ns~~~~~'Root:'~italic(F)~'= 0.050'^ns))+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 20,face="bold"),
        axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16, angle = -20),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20),
        plot.subtitle=element_text(face = "bold", colour = "black", size = 20))+
  facet_grid( .~env$Compartment)
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS AVD.pdf", p,
       height=8,width =9,limitsize = FALSE )
