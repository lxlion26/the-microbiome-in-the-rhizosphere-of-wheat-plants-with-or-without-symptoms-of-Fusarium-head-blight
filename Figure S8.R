###########################Fungal class barplot##################
setwd("C:/graduate/WH赤霉病/扩增子/WH Fungi")

rm(list = ls())
ID <- read.csv("ID.onlyFun.csv", head = T, row.names = 1)
ID$OTU<-rownames(ID)
data0 <- read.csv("da.rar.csv", head = T, row.names = 1)
env <- read.csv("env.csv", head = T)

ID$OTU.ID<-paste(rownames(ID),ID$phylum, ID$genus, sep="_")
sum(env$Sample != rownames(data0)) 
sum(rownames(ID) != colnames(data0))

str(ID$kingdom)

Flev<-ID[,"class"] 
fung.lev<-data.frame(aggregate(t(data0),by=list(Flev) , sum))
rownames(fung.lev)<-fung.lev[,1]; 
fung.lev<-fung.lev[,-1]
data1<-data.frame(t(fung.lev))

total<-apply(data1, 1, sum); 
fung.relabu<-data.frame(lapply(data1, function(x) {  x / total  }) )  # change the abundance data into relative abundance#

fung.RZ<-fung.relabu
env.RZ<-env

fung.RZ <- fung.RZ[,order(-colSums(fung.RZ))] ## 
lev<-interaction(env.RZ$Compartment, env.RZ$State,sep = ":") ## ccombining factors for Barplot profiling
fung.RZ.lev0<-aggregate(fung.RZ,by=list(lev) , mean) # generate the mean for each factor level
fung.RZ.lev <- fung.RZ.lev0[, !names(fung.RZ.lev0) %in% c("X_")]
fung.RZ1<-fung.RZ.lev[,c(1,2:11)] # the domiant OTUs
fung.RZ1 <- melt(fung.RZ1,id.vars = "Group.1")
fung.RZ1<-cSplit(fung.RZ1, "Group.1", ":")
names(fung.RZ1)<-c("Fungi","Relative_Abundance", "Compartment", "State")

fung.RZ2<-fung.RZ.lev[,c(1,12:ncol(fung.RZ.lev))] # the rare OTUs, combined as 'others'
fung.RZ2 <- melt(fung.RZ2,id.vars = "Group.1")
fung.RZ2<-cSplit(fung.RZ2, "Group.1", ":")
names(fung.RZ2)<-c("Fungi","Relative_Abundance", "Compartment", "State")
fung.RZ2$Fungi<-"Other"

fung.bind<-rbind(fung.RZ1,fung.RZ2) # combine the domiant and rare OTUs

col11<-c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#43A743","#97D1A0",  "grey","#ff00ff","#00ff00", "deepskyblue", "gold", "red", "navy", "darkgreen","maroon3", "black", "bisque", "grey")
p<-ggplot(fung.bind, aes(x = State, y = Relative_Abundance, fill=Fungi)) +
  geom_bar(stat='identity', position = "fill", width = 0.5)+
  facet_grid(~Compartment)+
  scale_fill_manual(values= col11)+
  theme_bw()+
  labs(x="",y = "Relative abundance")+
  guides(fill=guide_legend(title= "class"))+
  scale_y_continuous(labels = scales::percent)+
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 20,face="bold"),
        legend.title = element_text(colour="black", size=20, face="bold"),
        legend.text = element_text(colour="black", size=16),
        axis.text.y=element_text(colour="black",size=16),
        axis.text.x=element_text(colour="black",size=16,angle = -20),
        axis.title=element_text(colour="black",size=20,face="bold"))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS class barplot.pdf", p,
       height=8,width =9,limitsize = FALSE )

#########################Fungal genus barplot################
rm(list = ls())
ID <- read.csv("ID.onlyFun.csv", head = T, row.names = 1)
ID$OTU<-rownames(ID)
data0 <- read.csv("da.rar.csv", head = T, row.names = 1)
env <- read.csv("env.csv", head = T)

ID$OTU.ID<-paste(rownames(ID),ID$phylum, ID$genus, sep="_")
sum(env$Sample != rownames(data0)) 
sum(rownames(ID) != colnames(data0)) 

str(ID$kingdom)

Flev<-ID[,"genus"] 
fung.lev<-data.frame(aggregate(t(data0),by=list(Flev) , sum))
rownames(fung.lev)<-fung.lev[,1]; 
fung.lev<-fung.lev[,-1]
data1<-data.frame(t(fung.lev))

total<-apply(data1, 1, sum); 
fung.relabu<-data.frame(lapply(data1, function(x) {  x / total  }) )  # change the abundance data into relative abundance#

fung.RZ<-fung.relabu
env.RZ<-env

fung.RZ <- fung.RZ[,order(-colSums(fung.RZ))] ## 
lev<-interaction(env.RZ$Compartment, env.RZ$State,sep = ":") ## ccombining factors for Barplot profiling
fung.RZ.lev0<-aggregate(fung.RZ,by=list(lev) , mean) # generate the mean for each factor level
fung.RZ.lev <- fung.RZ.lev0[, !names(fung.RZ.lev0) %in% c("X_")]
fung.RZ1<-fung.RZ.lev[,c(1,2:11)] # the domiant OTUs
fung.RZ1 <- melt(fung.RZ1,id.vars = "Group.1")
fung.RZ1<-cSplit(fung.RZ1, "Group.1", ":")
names(fung.RZ1)<-c("Fungi","Relative_Abundance", "Compartment", "State")

fung.RZ2<-fung.RZ.lev[,c(1,12:ncol(fung.RZ.lev))] # the rare OTUs, combined as 'others'
fung.RZ2 <- melt(fung.RZ2,id.vars = "Group.1")
fung.RZ2<-cSplit(fung.RZ2, "Group.1", ":")
names(fung.RZ2)<-c("Fungi","Relative_Abundance", "Compartment", "State")
fung.RZ2$Fungi<-"Other"

fung.bind<-rbind(fung.RZ1,fung.RZ2) # combine the domiant and rare OTUs

col11<-c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#43A743","#97D1A0",  "grey","#ff00ff","#00ff00", "deepskyblue", "gold", "red", "navy", "darkgreen","maroon3", "black", "bisque", "grey")
p<-ggplot(fung.bind, aes(x = State, y = Relative_Abundance, fill=Fungi)) +
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
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS genus barplot.pdf", p,
       height=8,width =9,limitsize = FALSE )

#######################Fungi###########################
setwd("C:/graduate/WH赤霉病/扩增子/WH Fungi")
library(vegan)
library(ggplot2)
library(colorRamps)
library(ape)
library(splitstackshape)
library(reshape2)


rm(list = ls())
ID <- read.csv("ID.onlyFun.csv", head = T, row.names = 1)
ID$OTU<-rownames(ID)
data0 <- read.csv("da.rar.csv", head = T, row.names = 1)
env <- read.csv("env.csv", head = T)


ID$OTU.ID<-paste(rownames(ID),ID$phylum, ID$genus, sep="_")
sum(env$Sample != rownames(data0)) 
sum(rownames(ID) != colnames(data0))

Flev<-ID[,"class"] 
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

####################
##Sordariomycetes
#######################
shapiro.test(fung.R$Sordariomycetes)
shapiro.test(log(fung.R$Sordariomycetes))
shapiro.test(sqrt(fung.R$Sordariomycetes))
shapiro.test(fung.Z$Sordariomycetes)
shapiro.test(log(fung.Z$Sordariomycetes)) 
shapiro.test(sqrt(fung.Z$Sordariomycetes))

bartlett.test(fung.Z$Sordariomycetes ~ interaction(State),data=fung.Z)#??

anova_fung.Z_Sordariomycetes<- aov(fung.Z$Sordariomycetes ~ State,data = fung.Z)
summary(anova_fung.Z_Sordariomycetes)

kruskal.test(fung.R$Sordariomycetes ~ fung.R$State)#


p<-ggplot(fung.RZ, aes(x = State, y = Sordariomycetes, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Sordariomycetes",expression('Rhizosphere:'~italic(χ) ^2~'=  14.465***'~~~~~'Root:'~italic(F)~'= 4.091*'))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Sordariomycetes.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##Dothideomycetes
#######################

shapiro.test(fung.R$Dothideomycetes)
shapiro.test(log(fung.R$Dothideomycetes))
shapiro.test(sqrt(fung.R$Dothideomycetes))
shapiro.test(fung.Z$Dothideomycetes)
shapiro.test(log(fung.Z$Dothideomycetes)) 
shapiro.test(sqrt(fung.Z$Dothideomycetes))

bartlett.test(log(fung.R$Dothideomycetes) ~ interaction(State),data=fung.R)

anova_log_fung.R_Dothideomycetes<- aov(log(fung.R$Dothideomycetes) ~ State,data = fung.R)
summary(anova_log_fung.R_Dothideomycetes)

kruskal.test(fung.Z$Dothideomycetes ~ fung.Z$State)#


p<-ggplot(fung.RZ, aes(x = State, y = Dothideomycetes, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Dothideomycetes",expression('Rhizosphere:'~italic(F)~'=  11.340**'~~~~~'Root:'~italic(χ) ^2~'= 3.686'^ns))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Dothideomycetes.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##Tremellomycetes
#######################

shapiro.test(fung.R$Tremellomycetes)
shapiro.test(log(fung.R$Tremellomycetes))
shapiro.test(sqrt(fung.R$Tremellomycetes))
shapiro.test(fung.Z$Tremellomycetes)
shapiro.test(log(fung.Z$Tremellomycetes))
shapiro.test(sqrt(fung.Z$Tremellomycetes))


kruskal.test(fung.R$Tremellomycetes ~ fung.R$State)#
kruskal.test(fung.Z$Tremellomycetes ~ fung.Z$State)#

p<-ggplot(fung.RZ, aes(x = State, y = Tremellomycetes, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Tremellomycetes",expression('Rhizosphere:'~italic(χ) ^2~'=  5.924*'~~~~~'Root:'~italic(χ) ^2~'= 4.406*'))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Tremellomycetes.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##Leotiomycetes
#######################

shapiro.test(fung.R$Leotiomycetes)
shapiro.test(log(fung.R$Leotiomycetes))
shapiro.test(sqrt(fung.R$Leotiomycetes))
shapiro.test(fung.Z$Leotiomycetes)
shapiro.test(log(fung.Z$Leotiomycetes))
shapiro.test(sqrt(fung.Z$Leotiomycetes))

bartlett.test(log(fung.Z$Leotiomycetes) ~ interaction(State),data=fung.Z)#??

anova_log_fung.Z_Leotiomycetes<- aov(log(fung.Z$Leotiomycetes) ~ State,data = fung.Z)
summary(anova_log_fung.Z_Leotiomycetes)

kruskal.test(fung.R$Leotiomycetes ~ fung.R$State)#


p<-ggplot(fung.RZ, aes(x = State, y = Leotiomycetes, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Leotiomycetes",expression('Rhizosphere:'~italic(χ) ^2~'=  3.850*'~~~~~'Root:'~italic(F)~'= 2.942'^ns))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Leotiomycetes.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##Cystobasidiomycetes
#######################

shapiro.test(fung.R$Cystobasidiomycetes)
shapiro.test(log(fung.R$Cystobasidiomycetes))
shapiro.test(sqrt(fung.R$Cystobasidiomycetes))
shapiro.test(fung.Z$Cystobasidiomycetes)
shapiro.test(log(fung.Z$Cystobasidiomycetes))
shapiro.test(sqrt(fung.Z$Cystobasidiomycetes))


kruskal.test(fung.R$Cystobasidiomycetes ~ fung.R$State)#
kruskal.test(fung.Z$Cystobasidiomycetes ~ fung.Z$State)#


p<-ggplot(fung.RZ, aes(x = State, y = Cystobasidiomycetes, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Cystobasidiomycetes",expression('Rhizosphere:'~italic(χ) ^2~'=  4.445*'~~~~~'Root:'~italic(χ) ^2~'= 1.179'^ns))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Cystobasidiomycetes.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##Spizellomycetes
#######################

shapiro.test(fung.R$Spizellomycetes)#??
shapiro.test(log(fung.R$Spizellomycetes)) #??
shapiro.test(sqrt(fung.R$Spizellomycetes))# ??
shapiro.test(fung.Z$Spizellomycetes)#??
shapiro.test(log(fung.Z$Spizellomycetes)) #??
shapiro.test(sqrt(fung.Z$Spizellomycetes))#??


kruskal.test(fung.R$Spizellomycetes ~ fung.R$State)#
kruskal.test(fung.Z$Spizellomycetes ~ fung.Z$State)#


p<-ggplot(fung.RZ, aes(x = State, y = Spizellomycetes, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Spizellomycetes",expression('Rhizosphere:'~italic(χ) ^2~'=  6.235*'~~~~~'Root:'~italic(χ) ^2~'= 2.035'^ns))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Spizellomycetes.pdf", p,
       height=8,width =9,limitsize = FALSE )

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

####################
##Cladosporium
#######################

shapiro.test(fung.R$Cladosporium)
shapiro.test(log(fung.R$Cladosporium))
shapiro.test(sqrt(fung.R$Cladosporium))
shapiro.test(fung.Z$Cladosporium)
shapiro.test(log(fung.Z$Cladosporium))
shapiro.test(sqrt(fung.Z$Cladosporium))

bartlett.test(sqrt(fung.R$Cladosporium) ~ interaction(State),data=fung.R)


kruskal.test(fung.R$Cladosporium ~ fung.R$State)#
kruskal.test(fung.Z$Cladosporium ~ fung.Z$State)#

p<-ggplot(fung.RZ, aes(x = State, y = Cladosporium, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Cladosporium",expression('Rhizosphere:'~italic(χ) ^2~'=  6.119**'~~~~~'Root:'~italic(χ) ^2~'= 1.895'^ns))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Cladosporium.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##Acremonium
#######################

shapiro.test(fung.R$Acremonium)
shapiro.test(log(fung.R$Acremonium))
shapiro.test(sqrt(fung.R$Acremonium))
shapiro.test(fung.Z$Acremonium)
shapiro.test(log(fung.Z$Acremonium))
shapiro.test(sqrt(fung.Z$Acremonium))

bartlett.test(sqrt(fung.R$Acremonium) ~ interaction(State),data=fung.R)

anova_sqrt_fung.R_Acremonium<- aov(sqrt(fung.R$Acremonium) ~ State,data = fung.R)
summary(anova_sqrt_fung.R_Acremonium)

kruskal.test(fung.Z$Acremonium ~ fung.Z$State)#


p<-ggplot(fung.RZ, aes(x = State, y = Acremonium, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Acremonium",expression('Rhizosphere:'~italic(F)~'=  4.769*'~~~~~'Root:'~italic(χ) ^2~'= 4.185*'))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Acremonium.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##Oliveonia
#######################

shapiro.test(fung.R$Oliveonia)
shapiro.test(log(fung.R$Oliveonia))
shapiro.test(sqrt(fung.R$Oliveonia))
shapiro.test(fung.Z$Oliveonia)
shapiro.test(log(fung.Z$Oliveonia))
shapiro.test(sqrt(fung.Z$Oliveonia))


kruskal.test(fung.R$Oliveonia ~ fung.R$State)#
kruskal.test(fung.Z$Oliveonia ~ fung.Z$State)#


p<-ggplot(fung.RZ, aes(x = State, y = Oliveonia, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Oliveonia",expression('Rhizosphere:'~italic(χ) ^2~'=  8.651**'~~~~~'Root:'~italic(χ) ^2~'= 1.061'^ns))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Oliveonia.pdf", p,
       height=8,width =9,limitsize = FALSE )
####################
##Coprinellus
#######################

shapiro.test(fung.R$Coprinellus)
shapiro.test(log(fung.R$Coprinellus))
shapiro.test(sqrt(fung.R$Coprinellus))
shapiro.test(fung.Z$Coprinellus)
shapiro.test(log(fung.Z$Coprinellus))
shapiro.test(sqrt(fung.Z$Coprinellus))


kruskal.test(fung.R$Coprinellus ~ fung.R$State)#
kruskal.test(fung.Z$Coprinellus ~ fung.Z$State)#


p<-ggplot(fung.RZ, aes(x = State, y = Coprinellus, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Coprinellus",expression('Rhizosphere:'~italic(χ) ^2~'=  6.407*'~~~~~'Root:'~italic(χ) ^2~'= 5.122*'))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of Coprinellus.pdf", p,
       height=8,width =9,limitsize = FALSE )