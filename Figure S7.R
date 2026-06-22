######################Bacterial phylum barplot ##################
library(vegan)
library(ggplot2)
library(colorRamps)
library(ape)
library(splitstackshape)
library(reshape2)

setwd("C:/graduate/WH赤霉病/扩增子/WH Bacteria")
rm(list = ls())
ID <- read.csv("ID.onlyBac.csv", head = T, row.names = 1)
ID$OTU<-rownames(ID)
data0 <- read.csv("da.rar.csv", head = T, row.names = 1)
env <- read.csv("env.csv", head = T)

ID$OTU.ID<-paste(rownames(ID),ID$phylum, ID$genus, sep="_")
sum(env$Sample != rownames(data0))  
sum(rownames(ID) != colnames(data0))

str(ID$kingdom)

Blev<-ID[,"phylum"] 
bac.lev<-data.frame(aggregate(t(data0),by=list(Blev) , sum))
rownames(bac.lev)<-bac.lev[,1]; 
bac.lev<-bac.lev[,-1]
data1<-data.frame(t(bac.lev))

total<-apply(data1, 1, sum); 
bac.relabu<-data.frame(lapply(data1, function(x) {  x / total  }) )

bac.RZ<-bac.relabu
env.RZ<-env

bac.RZ <- bac.RZ[,order(-colSums(bac.RZ))] ## 
lev<-interaction(env.RZ$Compartment, env.RZ$State,sep = ":") ## ccombining factors for Barplot profiling
bac.RZ.lev<-aggregate(bac.RZ,by=list(lev) , mean) # generate the mean for each factor level
bac.RZ1<-bac.RZ.lev[,c(1,2:11)] # the domiant OTUs
bac.RZ1 <- melt(bac.RZ1,id.vars = "Group.1")
bac.RZ1<-cSplit(bac.RZ1, "Group.1", ":")

names(bac.RZ1)<-c("Bacteria","Relative_Abundance", "Compartment", "State")

bac.RZ2<-bac.RZ.lev[,c(1,12:ncol(bac.RZ.lev))] # the rare OTUs, combined as 'others'
bac.RZ2 <- melt(bac.RZ2,id.vars = "Group.1")
bac.RZ2<-cSplit(bac.RZ2, "Group.1", ":")
names(bac.RZ2)<-c("Bacteria","Relative_Abundance", "Compartment", "State")
bac.RZ2$Bacteria<-"Other"

bac.bind<-rbind(bac.RZ1,bac.RZ2) # combine the domiant and rare OTUs

col11<-c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#43A743","#97D1A0",  "grey","#ff00ff","#00ff00", "deepskyblue", "gold", "red", "navy", "darkgreen","maroon3", "black", "bisque", "grey")

p<-ggplot(bac.bind, aes(x = State, y = Relative_Abundance, fill=Bacteria)) +
  geom_bar(stat='identity', position = "fill", width = 0.5)+
  facet_grid(~Compartment)+
  scale_fill_manual(values= col11)+
  theme_bw()+
  labs(x="",y = "Relative abundance")+
  guides(fill=guide_legend(title= "phylum"))+
  scale_y_continuous(labels = scales::percent)+
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 20,face="bold"),
        legend.title = element_text(colour="black", size=20, face="bold"),
        legend.text = element_text(colour="black", size=16),
        axis.text.y=element_text(colour="black",size=16),
        axis.text.x=element_text(colour="black",size=16,angle = -20),
        axis.title=element_text(colour="black",size=20,face="bold"))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S phylum barplot.pdf", p,
       height=8,width =9,limitsize = FALSE )

#######################Bacterial genus barplot###############
rm(list = ls())
ID <- read.csv("ID.onlyBac.csv", head = T, row.names = 1)
ID$OTU<-rownames(ID)
data0 <- read.csv("da.rar.csv", head = T, row.names = 1)
env <- read.csv("env.csv", head = T)

ID$OTU.ID<-paste(rownames(ID),ID$phylum, ID$genus, sep="_")
sum(env$Sample != rownames(data0))  
sum(rownames(ID) != colnames(data0)) 

str(ID$kingdom)

Blev<-ID[,"genus"]
bac.lev<-data.frame(aggregate(t(data0),by=list(Blev) , sum))
rownames(bac.lev)<-bac.lev[,1]; 
bac.lev<-bac.lev[,-1]
data1<-data.frame(t(bac.lev))

total<-apply(data1, 1, sum); 
bac.relabu<-data.frame(lapply(data1, function(x) {  x / total  }) )

bac.RZ<-bac.relabu
env.RZ<-env

bac.RZ <- bac.RZ[,order(-colSums(bac.RZ))] ## 
lev<-interaction(env.RZ$Compartment, env.RZ$State,sep = ":") ## ccombining factors for Barplot profiling
bac.RZ.lev0<-aggregate(bac.RZ,by=list(lev) , mean) # generate the mean for each factor level
bac.RZ.lev <- bac.RZ.lev0[, !names(bac.RZ.lev0) %in% c("X_")]
bac.RZ1<-bac.RZ.lev[,c(1,2:11)] # the domiant OTUs
bac.RZ1 <- melt(bac.RZ1,id.vars = "Group.1")
bac.RZ1<-cSplit(bac.RZ1, "Group.1", ":")

names(bac.RZ1)<-c("Bacteria","Relative_Abundance", "Compartment", "State")

bac.RZ2<-bac.RZ.lev[,c(1,12:ncol(bac.RZ.lev))] # the rare OTUs, combined as 'others'
bac.RZ2 <- melt(bac.RZ2,id.vars = "Group.1")
bac.RZ2<-cSplit(bac.RZ2, "Group.1", ":")
names(bac.RZ2)<-c("Bacteria","Relative_Abundance", "Compartment", "State")
bac.RZ2$Bacteria<-"Other"

bac.bind<-rbind(bac.RZ1,bac.RZ2) # combine the domiant and rare OTUs

col11<-c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#43A743","#97D1A0",  "grey","#ff00ff","#00ff00", "deepskyblue", "gold", "red", "navy", "darkgreen","maroon3", "black", "bisque", "grey")

p<-ggplot(bac.bind, aes(x = State, y = Relative_Abundance, fill=Bacteria)) +
  geom_bar(stat='identity', position = "fill", width = 0.5)+
  facet_grid(~Compartment)+
  scale_fill_manual(values= col11)+
  theme_bw()+
  labs(x="",y = "Relative abundance")+
  guides(fill=guide_legend(title= "genus"))+
  scale_y_continuous(labels = scales::percent)+
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 20,face="bold"),
        legend.title = element_text(colour="black", size=20, face="bold"),
        legend.text = element_text(colour="black", size=16),
        axis.text.y=element_text(colour="black",size=16),
        axis.text.x=element_text(colour="black",size=16,angle = -20),
        axis.title=element_text(colour="black",size=20,face="bold"))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S genus barplot.pdf", p,
       height=8,width =12,limitsize = FALSE )

##################Bacteria#########################
library(vegan)
library(ggplot2)
library(colorRamps)
library(ape)
library(splitstackshape)
library(reshape2)

rm(list = ls())
ID <- read.csv("ID.onlyBac.csv", head = T, row.names = 1)
ID$OTU<-rownames(ID)
data0 <- read.csv("da.rar.csv", head = T, row.names = 1)
env <- read.csv("env.csv", head = T)


ID$OTU.ID<-paste(rownames(ID),ID$phylum, ID$genus, sep="_")
sum(env$Sample != rownames(data0))  
sum(rownames(ID) != colnames(data0))

Blev<-ID[,"phylum"] 
bac.lev<-data.frame(aggregate(t(data0),by=list(Blev) , sum))
rownames(bac.lev)<-bac.lev[,1]; 
bac.lev<-bac.lev[,-1]
data1<-data.frame(t(bac.lev))

total<-apply(data1, 1, sum); 
bac.relabu<-data.frame(lapply(data1, function(x) {  x / total  }) )

bac.RZ<-bac.relabu
env.RZ<-env
bac.RZ<-bac.RZ[,c(2,3,5,7,12,13,17,20,21,24)]
bac.RZ<-cbind(bac.RZ,env.RZ[,c(3,4)])
bac.R<-bac.RZ[bac.RZ$Compartment=="Rhizosphere",]
bac.Z<-bac.RZ[bac.RZ$Compartment=="Root",]

####################
##Actinobacteriota
#######################
##正态性检验，p>0.05表明符合正态性分布，符合用F检验，不符合用卡方检验。
shapiro.test(bac.R$Actinobacteriota)#×
shapiro.test(log(bac.R$Actinobacteriota)) #×
shapiro.test(sqrt(bac.R$Actinobacteriota))# ×
shapiro.test(bac.Z$Actinobacteriota)#×
shapiro.test(log(bac.Z$Actinobacteriota)) #×
shapiro.test(sqrt(bac.Z$Actinobacteriota))# 对
#检验方差是否齐次。p>0.05表明方差齐次。
bartlett.test(sqrt(bac.Z$Actinobacteriota) ~ interaction(State),data=bac.Z)#对
#方差齐次后，进行ANOVA方差分析
anova_sqrt_bac.Z_Actinobacteriota<- aov(sqrt(bac.Z$Actinobacteriota) ~ State,data = bac.Z)
summary(anova_sqrt_bac.Z_Actinobacteriota)
####不符合正态分布时，采用非参数检验(卡方检验)
kruskal.test(bac.R$Actinobacteriota ~ bac.R$State)#


p<-ggplot(bac.RZ, aes(x = State, y = Actinobacteriota, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Actinobacteriota",expression('Rhizosphere:'~italic(χ) ^2~'=  10.116**'~~~~~'Root:'~italic(F)~'= 0.683'^ns))+
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
        plot.subtitle=element_text(colour = "black", size = 20))+
  facet_grid( .~env$Compartment)
p
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Actinobacteriota.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##Verrucomicrobiota
#######################
##正态性检验，p>0.05表明符合正态性分布，符合用F检验，不符合用卡方检验。
shapiro.test(bac.R$Verrucomicrobiota)#×
shapiro.test(log(bac.R$Verrucomicrobiota)) #×
shapiro.test(sqrt(bac.R$Verrucomicrobiota))# ×
shapiro.test(bac.Z$Verrucomicrobiota)#×
shapiro.test(log(bac.Z$Verrucomicrobiota)) #×
shapiro.test(sqrt(bac.Z$Verrucomicrobiota))# 对

####不符合正态分布时，采用非参数检验(卡方检验)
kruskal.test(bac.R$Verrucomicrobiota ~ bac.R$State)#
kruskal.test(bac.Z$Verrucomicrobiota ~ bac.Z$State)

p<-ggplot(bac.RZ, aes(x = State, y = Verrucomicrobiota, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Verrucomicrobiota",expression('Rhizosphere:'~italic(χ) ^2~'=  4.325*'~~~~~'Root:'~italic(χ) ^2~'= 0.172'^ns))+
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
        plot.subtitle=element_text(colour = "black", size = 20))+
  facet_grid( .~env$Compartment)
p
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Verrucomicrobiota.pdf", p,
       height=8,width =9,limitsize = FALSE )


Blev<-ID[,"genus"] 
bac.lev<-data.frame(aggregate(t(data0),by=list(Blev) , sum))
rownames(bac.lev)<-bac.lev[,1]; 
bac.lev<-bac.lev[,-1]
data1<-data.frame(t(bac.lev))

total<-apply(data1, 1, sum); 
bac.relabu<-data.frame(lapply(data1, function(x) {  x / total  }) )

bac.RZ<-bac.relabu
env.RZ<-env
bac.RZ<-cbind(bac.RZ,env.RZ[,c(3,4)])
bac.R<-bac.RZ[bac.RZ$Compartment=="Rhizosphere",]
bac.Z<-bac.RZ[bac.RZ$Compartment=="Root",]

####################
##Pseudarthrobacter
#######################
##正态性检验，p>0.05表明符合正态性分布，符合用F检验，不符合用卡方检验。
shapiro.test(bac.R$Pseudarthrobacter)#×
shapiro.test(log(bac.R$Pseudarthrobacter)) #对
shapiro.test(sqrt(bac.R$Pseudarthrobacter))#对
shapiro.test(bac.Z$Pseudarthrobacter)#×
shapiro.test(log(bac.Z$Pseudarthrobacter)) #×
shapiro.test(sqrt(bac.Z$Pseudarthrobacter))#对
#检验方差是否齐次。p>0.05表明方差齐次。
bartlett.test(log(bac.R$Pseudarthrobacter) ~ interaction(State),data=bac.R)#对
bartlett.test(sqrt(bac.R$Pseudarthrobacter) ~ interaction(State),data=bac.R)#对
bartlett.test(sqrt(bac.Z$Pseudarthrobacter) ~ interaction(State),data=bac.Z)#对
#方差齐次后，进行ANOVA方差分析
anova_log_bac.R_Pseudarthrobacter<- aov(log(bac.R$Pseudarthrobacter) ~ State,data = bac.R)
summary(anova_log_bac.R_Pseudarthrobacter)
anova_sqrt_bac.R_Pseudarthrobacter<- aov(sqrt(bac.R$Pseudarthrobacter) ~ State,data = bac.R)
summary(anova_sqrt_bac.R_Pseudarthrobacter)
anova_sqrt_bac.Z_Pseudarthrobacter<- aov(sqrt(bac.Z$Pseudarthrobacter) ~ State,data = bac.Z)
summary(anova_sqrt_bac.Z_Pseudarthrobacter)


p<-ggplot(bac.RZ, aes(x = State, y = Pseudarthrobacter, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Pseudarthrobacter",expression('Rhizosphere:'~italic(F)~'=  18.760***'~~~~~'Root:'~italic(F)~'= 6.393*'))+
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
        plot.subtitle=element_text(colour = "black", size = 20))+
  facet_grid( .~env$Compartment)
p
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Pseudarthrobacter.pdf", p,
       height=8,width =9,limitsize = FALSE )


####################
##Pantoea
#######################
##正态性检验，p>0.05表明符合正态性分布，符合用F检验，不符合用卡方检验。
shapiro.test(bac.R$Pantoea)#×
shapiro.test(log(bac.R$Pantoea)) #×
shapiro.test(sqrt(bac.R$Pantoea))# ×
shapiro.test(bac.Z$Pantoea)#×
shapiro.test(log(bac.Z$Pantoea)) #对
shapiro.test(sqrt(bac.Z$Pantoea))#×
#检验方差是否齐次。p>0.05表明方差齐次。
bartlett.test(log(bac.Z$Pantoea) ~ interaction(State),data=bac.Z)#对
#方差齐次后，进行ANOVA方差分析
anova_log_bac.Z_Pantoea<- aov(log(bac.Z$Pantoea) ~ State,data = bac.Z)
summary(anova_log_bac.Z_Pantoea)
####不符合正态分布时，采用非参数检验(卡方检验)
kruskal.test(bac.R$Pantoea ~ bac.R$State)#


p<-ggplot(bac.RZ, aes(x = State, y = Pantoea, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Pantoea",expression('Rhizosphere:'~italic(χ) ^2~'=  11.511***'~~~~~'Root:'~italic(F)~'= 2.472'^ns))+
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
        plot.subtitle=element_text(colour = "black", size = 20))+
  facet_grid( .~env$Compartment)
p
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Pantoea.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##Bacillus
#######################
##正态性检验，p>0.05表明符合正态性分布，符合用F检验，不符合用卡方检验。
shapiro.test(bac.R$Bacillus)#×
shapiro.test(log(bac.R$Bacillus)) #×
shapiro.test(sqrt(bac.R$Bacillus))# ×
shapiro.test(bac.Z$Bacillus)#×
shapiro.test(log(bac.Z$Bacillus)) #×
shapiro.test(sqrt(bac.Z$Bacillus))#×

####不符合正态分布时，采用非参数检验(卡方检验)
kruskal.test(bac.R$Bacillus ~ bac.R$State)#
kruskal.test(bac.Z$Bacillus ~ bac.Z$State)#


p<-ggplot(bac.RZ, aes(x = State, y = Bacillus, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Bacillus",expression('Rhizosphere:'~italic(χ) ^2~'=  4.735*'~~~~~'Root:'~italic(χ) ^2~'= 0.028'^ns))+
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
        plot.subtitle=element_text(colour = "black", size = 20))+
  facet_grid( .~env$Compartment)
p
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Bacillus.pdf", p,
       height=8,width =9,limitsize = FALSE )