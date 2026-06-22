setwd("C:/graduate/WH赤霉病/宏基因组/kraken")
library(vegan)
library(ggplot2)
library(colorRamps)
library(ape)
library(splitstackshape)
library(reshape2)

rm(list = ls())
ID <- read.csv("Fun.ID.csv", head = T, row.names = 1)
ID$OTU<-rownames(ID)
data0 <- read.csv("Fun.da.rar.csv", head = T, row.names = 1)
env <- read.csv("env.csv", head = T)


sum(env$Sample != rownames(data0))  
sum(ID$OTU.ID != colnames(data0))

Flev<-ID[,"species"] 
fun.lev<-data.frame(aggregate(t(data0),by=list(Flev) , sum))
rownames(fun.lev)<-fun.lev[,1]; 
fun.lev<-fun.lev[,-1]
data1<-data.frame(t(fun.lev))

total<-apply(data1, 1, sum); 
fun.relabu<-data.frame(lapply(data1, function(x) {  x / total  }) )
fun.R<-cbind(fun.relabu,env[,c(3,4)])
####################
##Fusarium falciforme
#######################
##正态性检验，p>0.05表明符合正态性分布，符合用F检验，不符合用卡方检验。
shapiro.test(fun.R$falciforme)#×
shapiro.test(log(fun.R$falciforme)) #×
shapiro.test(sqrt(fun.R$falciforme))# ×

####不符合正态分布时，采用非参数检验(卡方检验)
kruskal.test(fun.R$falciforme~ env$State)#


p<-ggplot(fun.R, aes(x = State, y = falciforme, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Fusarium falciforme",expression('Rhizosphere:'~italic(χ) ^2~'=  1.488'^ns~~~df~'=1'~~~italic(P)~'= 0.223'))+
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
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/RA of Fusarium falciforme.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##Fusarium keratoplasticum
#######################
##正态性检验，p>0.05表明符合正态性分布，符合用F检验，不符合用卡方检验。
shapiro.test(fun.R$keratoplasticum)#×
shapiro.test(log(fun.R$keratoplasticum)) #yes
shapiro.test(sqrt(fun.R$keratoplasticum))# ×

#检验方差是否齐次。p>0.05表明方差齐次。
bartlett.test(log(fun.R$keratoplasticum) ~ interaction(State),data=fun.R)#对
#方差齐次后，进行ANOVA方差分析
anova_log_fun.R_keratoplasticum<- aov(log(fun.R$keratoplasticum) ~ State,data = fun.R)
summary(anova_log_fun.R_keratoplasticum)

p<-ggplot(fun.R, aes(x = State, y = keratoplasticum, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Fusarium keratoplasticum",expression('Rhizosphere:'~italic(F)~'=  0.326'^ns~~~df~'=1'~~~italic(P)~'= 0.570'))+
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
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/RA of Fusarium keratoplasticum.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##Fusarium pseudograminearum
#######################
##正态性检验，p>0.05表明符合正态性分布，符合用F检验，不符合用卡方检验。
shapiro.test(fun.R$pseudograminearum)#×
shapiro.test(log(fun.R$pseudograminearum)) #×
shapiro.test(sqrt(fun.R$pseudograminearum))# ×

####不符合正态分布时，采用非参数检验(卡方检验)
kruskal.test(fun.R$pseudograminearum~ env$State)#

p<-ggplot(fun.R, aes(x = State, y = pseudograminearum, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Fusarium pseudograminearum",expression('Rhizosphere:'~italic(χ) ^2~'=  0.825'^ns~~~df~'=1'~~~italic(P)~'= 0.364'))+
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
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/RA of Fusarium pseudograminearum.pdf", p,
       height=8,width =9,limitsize = FALSE )


####################
##Fusarium poae
#######################
##正态性检验，p>0.05表明符合正态性分布，符合用F检验，不符合用卡方检验。
shapiro.test(fun.R$poae)#yes
shapiro.test(log(fun.R$poae)) #yes
shapiro.test(sqrt(fun.R$poae))# yes

#检验方差是否齐次。p>0.05表明方差齐次。
bartlett.test(fun.R$poae ~ interaction(State),data=fun.R)#yes
bartlett.test(log(fun.R$poae) ~ interaction(State),data=fun.R)#yes
bartlett.test(sqrt(fun.R$poae) ~ interaction(State),data=fun.R)#yes

#方差齐次后，进行ANOVA方差分析
anova_fun.R_poae<- aov(fun.R$poae ~ State,data = fun.R)
summary(anova_fun.R_poae)
anova_log_fun.R_poae<- aov(log(fun.R$poae) ~ State,data = fun.R)
summary(anova_log_fun.R_poae)
anova_sqrt_fun.R_poae<- aov(sqrt(fun.R$poae) ~ State,data = fun.R)
summary(anova_sqrt_fun.R_poae)


p<-ggplot(fun.R, aes(x = State, y = poae, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Fusarium poae",expression('Rhizosphere:'~italic(F)~'=  0.018'^ns~~~df~'=1'~~~italic(P)~'= 0.894'))+
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
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/RA of Fusarium poae.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##Fusarium graminearum
#######################
##正态性检验，p>0.05表明符合正态性分布，符合用F检验，不符合用卡方检验。
shapiro.test(fun.R$graminearum)#×
shapiro.test(log(fun.R$graminearum)) #×
shapiro.test(sqrt(fun.R$graminearum))# ×

####不符合正态分布时，采用非参数检验(卡方检验)
kruskal.test(fun.R$graminearum~ env$State)#


p<-ggplot(fun.R, aes(x = State, y = graminearum, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Fusarium graminearum",expression('Rhizosphere:'~italic(χ) ^2~'=  5.058*'~~~df~'=1'~~~italic(P)~'= 0.025'))+
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
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/RA of Fusarium graminearum.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##Fusarium oxysporum
#######################
##正态性检验，p>0.05表明符合正态性分布，符合用F检验，不符合用卡方检验。
shapiro.test(fun.R$oxysporum)#×
shapiro.test(log(fun.R$oxysporum)) #×
shapiro.test(sqrt(fun.R$oxysporum))# ×

####不符合正态分布时，采用非参数检验(卡方检验)
kruskal.test(fun.R$oxysporum~ env$State)#


p<-ggplot(fun.R, aes(x = State, y = oxysporum, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Fusarium oxysporum",expression('Rhizosphere:'~italic(χ) ^2~'=  6.163**'~~~df~'=1'~~~italic(P)~'= 0.013'))+
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
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/RA of Fusarium oxysporum.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##Fusarium fujikuroi
#######################
##正态性检验，p>0.05表明符合正态性分布，符合用F检验，不符合用卡方检验。
shapiro.test(fun.R$fujikuroi)#×
shapiro.test(log(fun.R$fujikuroi)) #×
shapiro.test(sqrt(fun.R$fujikuroi))# ×

####不符合正态分布时，采用非参数检验(卡方检验)
kruskal.test(fun.R$fujikuroi~ env$State)#


p<-ggplot(fun.R, aes(x = State, y = fujikuroi, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Fusarium fujikuroi",expression('Rhizosphere:'~italic(χ) ^2~'=  1.844'^ns~~~df~'=1'~~~italic(P)~'= 0.174'))+
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
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/RA of Fusarium fujikuroi.pdf", p,
       height=8,width =9,limitsize = FALSE )



####################
##Fusarium verticillioides
#######################
##正态性检验，p>0.05表明符合正态性分布，符合用F检验，不符合用卡方检验。
shapiro.test(fun.R$verticillioides)#×
shapiro.test(log(fun.R$verticillioides)) #yes
shapiro.test(sqrt(fun.R$verticillioides))# ×

#检验方差是否齐次。p>0.05表明方差齐次。
bartlett.test(log(fun.R$verticillioides) ~ interaction(State),data=fun.R)#yes

#方差齐次后，进行ANOVA方差分析
anova_log_fun.R_verticillioides<- aov(log(fun.R$verticillioides) ~ State,data = fun.R)
summary(anova_log_fun.R_verticillioides)

####不符合正态分布时，采用非参数检验(卡方检验)
kruskal.test(fun.R$verticillioides~ env$State)#

p<-ggplot(fun.R, aes(x = State, y = verticillioides, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Fusarium verticillioides",expression('Rhizosphere:'~italic(F)~'=  6.108*'~~~df~'=1'~~~italic(P)~'= 0.017'))+
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
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/RA of Fusarium verticillioides.pdf", p,
       height=8,width =9,limitsize = FALSE )

####################
##Fusarium musae
#######################
##正态性检验，p>0.05表明符合正态性分布，符合用F检验，不符合用卡方检验。
shapiro.test(fun.R$musae)#×
shapiro.test(log(fun.R$musae)) # ×
shapiro.test(sqrt(fun.R$musae))# ×


####不符合正态分布时，采用非参数检验(卡方检验)
kruskal.test(fun.R$musae~ env$State)#

p<-ggplot(fun.R, aes(x = State, y = musae, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Fusarium musae",expression('Rhizosphere:'~italic(χ) ^2~'=  5.176*'~~~df~'=1'~~~italic(P)~'= 0.023'))+
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
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/RA of Fusarium musae.pdf", p,
       height=8,width =9,limitsize = FALSE )


####################
##Fusarium 
#######################
rm(list = ls())
ID <- read.csv("Fun.ID.csv", head = T, row.names = 1)
ID$OTU<-rownames(ID)
data0 <- read.csv("Fun.da.rar.csv", head = T, row.names = 1)
env <- read.csv("env.csv", head = T)


sum(env$Sample != rownames(data0))  
sum(ID$OTU.ID != colnames(data0))

Flev<-ID[,"genus"] 
fun.lev<-data.frame(aggregate(t(data0),by=list(Flev) , sum))
rownames(fun.lev)<-fun.lev[,1]; 
fun.lev<-fun.lev[,-1]
data1<-data.frame(t(fun.lev))

total<-apply(data1, 1, sum); 
fun.relabu<-data.frame(lapply(data1, function(x) {  x / total  }) )
fun.R<-cbind(fun.relabu,env[,c(3,4)])

##正态性检验，p>0.05表明符合正态性分布，符合用F检验，不符合用卡方检验。
shapiro.test(fun.R$Fusarium)#×
shapiro.test(log(fun.R$Fusarium)) # yes
shapiro.test(sqrt(fun.R$Fusarium))# ×

#检验方差是否齐次。p>0.05表明方差齐次。
bartlett.test(log(fun.R$Fusarium) ~ interaction(State),data=fun.R)#yes

#方差齐次后，进行ANOVA方差分析
anova_log_fun.R_Fusarium<- aov(log(fun.R$Fusarium) ~ State,data = fun.R)
summary(anova_log_fun.R_Fusarium)

p<-ggplot(fun.R, aes(x = State, y = Fusarium, fill = State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.2,size=1.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  labs(title="",x = '', y = 'Relative abundance')+
  ggtitle("Fusarium",expression('Rhizosphere:'~italic(F)~'=  1.565'^ns~~~df~'=1'~~~italic(P)~'= 0.216'))+
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
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/RA of Fusarium.pdf", p,
       height=8,width =9,limitsize = FALSE )
