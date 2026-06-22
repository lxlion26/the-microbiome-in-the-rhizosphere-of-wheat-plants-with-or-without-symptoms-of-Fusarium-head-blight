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
sum(env$Sample != rownames(data0))  # 0 ˵??env??Otu table ??sample????˳????һ?µ?
sum(rownames(ID) != colnames(data0)) # 0 ˵??ID??Otu table ??OTU????˳????һ?µ?
#rownames(data0) <- ID0$ID.all


#############class###############
Flev<-ID[,"OTU"] 
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
##OTU_Fun5
#######################
##??̬?Լ??飬p>0.05??????????̬?Էֲ?????????F???飬???????ÿ??????顣
shapiro.test(fung.R$OTU5)#no
shapiro.test(log(fung.R$OTU5)) #no
shapiro.test(sqrt(fung.R$OTU5))# no
shapiro.test(fung.Z$OTU5)#no
shapiro.test(log(fung.Z$OTU5)) #yes
shapiro.test(sqrt(fung.Z$OTU5))#no
#???鷽???Ƿ????Ρ?p>0.05???????????Ρ?
bartlett.test(log(fung.Z$OTU5) ~ interaction(State),data=fung.Z)#yes
#???????κ󣬽???ANOVA????????
anova_fung.Z_log_OTU5<- aov(log(fung.Z$OTU5) ~ State,data = fung.Z)
summary(anova_fung.Z_log_OTU5)
####????????̬?ֲ?ʱ?????÷ǲ???????(????????)
kruskal.test(fung.R$OTU5 ~ fung.R$State)#


p<-ggplot(fung.RZ, aes(x = State, y = OTU5, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("OTU_Fun5",expression('Rhizosphere:'~italic(χ) ^2~'=  1.202'^ns~~~~~'Root:'~italic(F)~'= 0.360'^ns))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of OTU_Fun5.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##OTU_Fun15
#######################
##??̬?Լ??飬p>0.05??????????̬?Էֲ?????????F???飬???????ÿ??????顣
shapiro.test(fung.R$OTU15)#no
shapiro.test(log(fung.R$OTU15)) #no
shapiro.test(sqrt(fung.R$OTU15))# no
shapiro.test(fung.Z$OTU15)#no
shapiro.test(log(fung.Z$OTU15)) #no
shapiro.test(sqrt(fung.Z$OTU15))#no

####????????̬?ֲ?ʱ?????÷ǲ???????(????????)
kruskal.test(fung.R$OTU15 ~ fung.R$State)#
kruskal.test(fung.Z$OTU15 ~ fung.Z$State)


p<-ggplot(fung.RZ, aes(x = State, y = OTU15, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("OTU_Fun15",expression('Rhizosphere:'~italic(χ) ^2~'=  0.126'^ns~~~~~'Root:'~italic(χ) ^2~'= 0.138'^ns))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of OTU_Fun15.pdf", p,
       height=8,width =9,limitsize = FALSE )


####################
##OTU_Fun21
#######################
##??̬?Լ??飬p>0.05??????????̬?Էֲ?????????F???飬???????ÿ??????顣
shapiro.test(fung.R$OTU21)#no
shapiro.test(log(fung.R$OTU21)) #no
shapiro.test(sqrt(fung.R$OTU21))# no
shapiro.test(fung.Z$OTU21)#no
shapiro.test(log(fung.Z$OTU21)) #no
shapiro.test(sqrt(fung.Z$OTU21))#no

####????????̬?ֲ?ʱ?????÷ǲ???????(????????)
kruskal.test(fung.R$OTU21 ~ fung.R$State)#
kruskal.test(fung.Z$OTU21 ~ fung.Z$State)


p<-ggplot(fung.RZ, aes(x = State, y = OTU21, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("OTU_Fun21",expression('Rhizosphere:'~italic(χ) ^2~'=  0.204'^ns~~~~~'Root:'~italic(χ) ^2~'= 1.514'^ns))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of OTU_Fun21.pdf", p,
       height=8,width =9,limitsize = FALSE )


####################
##OTU_Fun29
#######################
##??̬?Լ??飬p>0.05??????????̬?Էֲ?????????F???飬???????ÿ??????顣
shapiro.test(fung.R$OTU29)#no
shapiro.test(log(fung.R$OTU29)) #yes
shapiro.test(sqrt(fung.R$OTU29))# no
shapiro.test(fung.Z$OTU29)#no
shapiro.test(log(fung.Z$OTU29)) #no
shapiro.test(sqrt(fung.Z$OTU29))#no
#???鷽???Ƿ????Ρ?p>0.05???????????Ρ?
bartlett.test(log(fung.R$OTU29) ~ interaction(State),data=fung.R)#yes
#???????κ󣬽???ANOVA????????
anova_fung.R_log_OTU29<- aov(log(fung.R$OTU29) ~ State,data = fung.R)
summary(anova_fung.R_log_OTU29)
####????????̬?ֲ?ʱ?????÷ǲ???????(????????)
kruskal.test(fung.Z$OTU29 ~ fung.Z$State)#

p<-ggplot(fung.RZ, aes(x = State, y = OTU21, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("OTU_Fun29",expression('Rhizosphere:'~italic(F)~'=  0.053'^ns~~~~~'Root:'~italic(χ) ^2~'= 0.840'^ns))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of OTU_Fun29.pdf", p,
       height=8,width =9,limitsize = FALSE )


####################
##OTU_Fun137
#######################
##??̬?Լ??飬p>0.05??????????̬?Էֲ?????????F???飬???????ÿ??????顣
shapiro.test(fung.R$OTU137)#no
shapiro.test(log(fung.R$OTU137)) #no
shapiro.test(sqrt(fung.R$OTU137))# no
shapiro.test(fung.Z$OTU137)#no
shapiro.test(log(fung.Z$OTU137)) #no
shapiro.test(sqrt(fung.Z$OTU137))#no

####????????̬?ֲ?ʱ?????÷ǲ???????(????????)
kruskal.test(fung.R$OTU137 ~ fung.R$State)#
kruskal.test(fung.Z$OTU137 ~ fung.Z$State)

p<-ggplot(fung.RZ, aes(x = State, y = OTU137, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("OTU_Fun137",expression('Rhizosphere:'~italic(χ) ^2~'=  0.020'^ns~~~~~'Root:'~italic(χ) ^2~'= 2.849'^ns))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of OTU_Fun137.pdf", p,
       height=8,width =9,limitsize = FALSE )


####################
##OTU_Fun535
#######################
##??̬?Լ??飬p>0.05??????????̬?Էֲ?????????F???飬???????ÿ??????顣
shapiro.test(fung.R$OTU535)#no
shapiro.test(log(fung.R$OTU535)) #no
shapiro.test(sqrt(fung.R$OTU535))# no
shapiro.test(fung.Z$OTU535)#no
shapiro.test(log(fung.Z$OTU535)) #no
shapiro.test(sqrt(fung.Z$OTU535))#no

####????????̬?ֲ?ʱ?????÷ǲ???????(????????)
kruskal.test(fung.R$OTU535 ~ fung.R$State)#
kruskal.test(fung.Z$OTU535 ~ fung.Z$State)

p<-ggplot(fung.RZ, aes(x = State, y = OTU535, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("OTU_Fun535",expression('Rhizosphere:'~italic(χ) ^2~'= 4.091*'~~~~~'Root:'~italic(χ) ^2~'= 0.001'^ns))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of OTU_Fun535.pdf", p,
       height=8,width =9,limitsize = FALSE )


####################
##OTU_Fun970
#######################
##??̬?Լ??飬p>0.05??????????̬?Էֲ?????????F???飬???????ÿ??????顣
shapiro.test(fung.R$OTU970)#no
shapiro.test(log(fung.R$OTU970)) #no
shapiro.test(sqrt(fung.R$OTU970))# no
shapiro.test(fung.Z$OTU970)#no
shapiro.test(log(fung.Z$OTU970)) #no
shapiro.test(sqrt(fung.Z$OTU970))#no

####????????̬?ֲ?ʱ?????÷ǲ???????(????????)
kruskal.test(fung.R$OTU970 ~ fung.R$State)#
kruskal.test(fung.Z$OTU970 ~ fung.Z$State)

p<-ggplot(fung.R, aes(x = State, y = OTU970, fill = State)) +
  geom_jitter(aes(color = State),position =position_jitter(0.2),alpha=0.8,size=2) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("OTU_Fun970",expression('Rhizosphere:'~italic(χ) ^2~'= 1.000'^ns))+
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
        plot.subtitle=element_text(colour = "black", size = 20))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of OTU_Fun970.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##OTU_Fun1057
#######################
##??̬?Լ??飬p>0.05??????????̬?Էֲ?????????F???飬???????ÿ??????顣
shapiro.test(fung.R$OTU1057)#no
shapiro.test(log(fung.R$OTU1057)) #no
shapiro.test(sqrt(fung.R$OTU1057))# no
shapiro.test(fung.Z$OTU1057)#no
shapiro.test(log(fung.Z$OTU1057)) #no
shapiro.test(sqrt(fung.Z$OTU1057))#no

####????????̬?ֲ?ʱ?????÷ǲ???????(????????)
kruskal.test(fung.R$OTU1057 ~ fung.R$State)#
kruskal.test(fung.Z$OTU1057 ~ fung.Z$State)

p<-ggplot(fung.RZ, aes(x = State, y = OTU1057, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("OTU_Fun1057",expression('Rhizosphere:'~italic(χ) ^2~'= 3.328'^ns~~~~~'Root:'~italic(χ) ^2~'= 2.424'^ns))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of OTU_Fun1057.pdf", p,
       height=8,width =9,limitsize = FALSE )


####################
##OTU_Fun1373
#######################
##??̬?Լ??飬p>0.05??????????̬?Էֲ?????????F???飬???????ÿ??????顣
shapiro.test(fung.R$OTU1373)#no
shapiro.test(log(fung.R$OTU1373)) #no
shapiro.test(sqrt(fung.R$OTU1373))# no
shapiro.test(fung.Z$OTU1373)#no
shapiro.test(log(fung.Z$OTU1373)) #no
shapiro.test(sqrt(fung.Z$OTU1373))#no


####????????̬?ֲ?ʱ?????÷ǲ???????(????????)
kruskal.test(fung.R$OTU1373 ~ fung.R$State)#
kruskal.test(fung.Z$OTU1373 ~ fung.Z$State)


p<-ggplot(fung.RZ, aes(x = State, y = OTU1373, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("OTU_Fun1373",expression('Rhizosphere:'~italic(χ) ^2~'= 2.022'^ns~~~~~'Root:'~italic(χ) ^2~'= 1.483'^ns))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of OTU_Fun1373.pdf", p,
       height=8,width =9,limitsize = FALSE )


####################
##OTU_Fun1438
#######################
##??̬?Լ??飬p>0.05??????????̬?Էֲ?????????F???飬???????ÿ??????顣
shapiro.test(fung.R$OTU1438)#no
shapiro.test(log(fung.R$OTU1438)) #no
shapiro.test(sqrt(fung.R$OTU1438))# no
shapiro.test(fung.Z$OTU1438)#no
shapiro.test(log(fung.Z$OTU1438)) #no
shapiro.test(sqrt(fung.Z$OTU1438))#no

####????????̬?ֲ?ʱ?????÷ǲ???????(????????)
kruskal.test(fung.R$OTU1438 ~ fung.R$State)#
kruskal.test(fung.Z$OTU1438 ~ fung.Z$State)

p<-ggplot(fung.RZ, aes(x = State, y = OTU1438, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("OTU_Fun1438",expression('Rhizosphere:'~italic(χ) ^2~'= 2.065'^ns~~~~~'Root:'~italic(χ) ^2~'= 7.412**'))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of OTU_Fun1438.pdf", p,
       height=8,width =9,limitsize = FALSE )



####################
##OTU_Fun1782
#######################
##??̬?Լ??飬p>0.05??????????̬?Էֲ?????????F???飬???????ÿ??????顣
shapiro.test(fung.R$OTU1782)#no
shapiro.test(log(fung.R$OTU1782)) #no
shapiro.test(sqrt(fung.R$OTU1782))# no
shapiro.test(fung.Z$OTU1782)#no
shapiro.test(log(fung.Z$OTU1782)) #no
shapiro.test(sqrt(fung.Z$OTU1782))#no


####????????̬?ֲ?ʱ?????÷ǲ???????(????????)
kruskal.test(fung.R$OTU1782 ~ fung.R$State)#
kruskal.test(fung.Z$OTU1782 ~ fung.Z$State)

p<-ggplot(fung.Z, aes(x = State, y = OTU1782, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("OTU_Fun1782",expression('Root:'~italic(χ) ^2~'= 0.178'^ns))+
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
        plot.subtitle=element_text(colour = "black", size = 20))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of OTU_Fun1782.pdf", p,
       height=8,width =9,limitsize = FALSE )


####################
##OTU_Fun1786
#######################
##??̬?Լ??飬p>0.05??????????̬?Էֲ?????????F???飬???????ÿ??????顣
shapiro.test(fung.R$OTU1786)#no
shapiro.test(log(fung.R$OTU1786)) #no
shapiro.test(sqrt(fung.R$OTU1786))# no
shapiro.test(fung.Z$OTU1786)#no
shapiro.test(log(fung.Z$OTU1786)) #no
shapiro.test(sqrt(fung.Z$OTU1786))#no

####????????̬?ֲ?ʱ?????÷ǲ???????(????????)
kruskal.test(fung.R$OTU1786 ~ fung.R$State)#
kruskal.test(fung.Z$OTU1786 ~ fung.Z$State)


p<-ggplot(fung.R, aes(x = State, y = OTU1786, fill = State)) +
  geom_jitter(aes(color = State),position =position_jitter(0.2),alpha=0.8,size=2) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("OTU_Fun1786",expression('Rhizosphere:'~italic(χ) ^2~'= 1.000'^ns))+
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
        plot.subtitle=element_text(colour = "black", size = 20))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of OTU_Fun1786.pdf", p,
       height=8,width =9,limitsize = FALSE )


####################
##OTU_Fun1903
#######################
##??̬?Լ??飬p>0.05??????????̬?Էֲ?????????F???飬???????ÿ??????顣
shapiro.test(fung.R$OTU1903)#no
shapiro.test(log(fung.R$OTU1903)) #no
shapiro.test(sqrt(fung.R$OTU1903))# no
shapiro.test(fung.Z$OTU1903)#no
shapiro.test(log(fung.Z$OTU1903)) #no
shapiro.test(sqrt(fung.Z$OTU1903))#no

####????????̬?ֲ?ʱ?????÷ǲ???????(????????)
kruskal.test(fung.R$OTU1903 ~ fung.R$State)#
kruskal.test(fung.Z$OTU1903 ~ fung.Z$State)

p<-ggplot(fung.RZ, aes(x = State, y = OTU1903, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("OTU_Fun1903",expression('Rhizosphere:'~italic(χ) ^2~'= 1.390'^ns~~~~~'Root:'~italic(χ) ^2~'= 4.252*'))+
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
ggsave("C:/graduate/WH赤霉病/pdf picture/RA of OTU_Fun1903.pdf", p,
       height=8,width =9,limitsize = FALSE )
