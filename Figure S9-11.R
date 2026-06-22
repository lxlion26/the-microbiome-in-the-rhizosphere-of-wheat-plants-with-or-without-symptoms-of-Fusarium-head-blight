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

################Bacterial LefSe##########################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
library(coin) # for the kruskal_test and wilcox_test
library(MicrobiotaProcess)
library(phyloseq)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(ggstar)
library(forcats)

setwd("C:/graduate/WH赤霉病/扩增子/WH Bacteria")

rm(list=ls())
otumat <- read.csv('da.rar.csv', row.names = 1)
otumat<-as.data.frame(t(otumat))
otumat = data.matrix(otumat)

taxmat <- read.csv('LEfSe细菌ID.csv', row.names = 1)
taxmat = as.matrix(taxmat)

sampledata <- read.csv('细菌LEfSe分组信息.csv', row.names = 1)
sampledata = sample_data(data.frame(sampledata))


OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(taxmat)
OTU
TAX
sampledata
#加入phyloseq体系
physeq = phyloseq(OTU, TAX,sampledata)
physeq


set.seed(1024)#由于LDA效应大小是通过随机重采样计算的，因此应设置随机种子以实现结果的重复性
deres <- diff_analysis(obj = physeq, 
                       classgroup = "Group",#分组
                       firstcomfun = "kruskal_test",
                       padjust = "fdr",#p值校正方法
                       filtermod = "pvalue",#以pvalue列过滤
                       firstalpha = 0.05,
                       strictmod = TRUE,#是否进行一对一事后检验
                       secondcomfun = "wilcox_test",
                       secondalpha = 0.01,
                       mlfun = "lda",#线性判别分析，可选随机森林
                       ldascore=3#线性判别分数
)
deres


p<-ggdiffclade(obj=deres,
            layout="circular",#布局类型
            alpha=0.3, #树分支的背景透明度
            linewd=0.2, #树中连线粗细
            skpointsize=0.8, #树骨架中国点的大小
            taxlevel=2, #要展示的树分支水平2（门）及3以下（纲目科属种）
            cladetext=4,#文本大小
            setColors=F,#自定义颜色
            removeUnknown=T,#不删分支，但移除分类中有un_的物种注释
            reduce=T)+ # 移除分类不明的物种分支和其注释，为T时removeUnknown参数无效
  scale_fill_manual(values=c("#8b66b8","#51c4c2")) +
  guides(color = guide_legend(keywidth = 0.1, keyheight = 0.6,order = 3,ncol=1)) +
  theme(panel.background=element_rect(fill=NA),
        legend.position="right", 
        plot.margin=margin(0,0,0,0),
        legend.spacing.y=unit(0.1, "cm"), 
        legend.title=element_text(size=20),
        legend.text=element_text(size=16), 
        legend.box.spacing=unit(0.1,"cm"))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S LefSe(1).pdf", p,
       height=8,width =12,limitsize = FALSE )

# visualization of different results by ggdiffbox
diffbox <- ggdiffbox(obj=deres, box_notch=FALSE, 
                     colorlist=c("#8b66b8","#c63596","#51c4c2","#2AB34A"), 
                     l_xlabtext="Relative abundance")
diffbox
ggsave("C:/graduate/WH赤霉病/pdf picture/16S LefSe biomarker(1).pdf", diffbox,
       height=8,width =9,limitsize = FALSE )

effectsize <- ggeffectsize(obj=deres, lineheight=0.1,linewidth=0.3,pointsize=3) + 
  scale_color_manual(values=c("#8b66b8","#51c4c2")) 
effectsize
ggsave("C:/graduate/WH赤霉病/pdf picture/16S LefSe biomarker(2).pdf", effectsize,
       height=8,width =9,limitsize = FALSE )

ggdifftaxbar(obj=deres, xtextsize=1.5,figwidth = 10, ylabel = "Relative abundance",output="each_biomarker_barplot",coloslist=c("#8b66b8","#c63596","#51c4c2","#2AB34A"))


sss <- as.MPSE(physeq)
sss %<>% mp_rrarefy() 
sss %<>% mp_cal_abundance(.abundance = RareAbundance,force=T) %>% mp_cal_abundance(.abundance=RareAbundance,.group=Group,force = T)
sss %<>% mp_diff_analysis(
  .abundance = RelRareAbundanceBySample,
  .group = Group,
  p.adjust="fdr",
  filter.p="pvalue",
  first.test.alpha = 0.05,
  second.test.alpha = 0.01,
  ldascore = 3
)
taxa.tree <- sss %>% mp_extract_tree(type="taxatree")
taxa.tree %>% select(label, nodeClass, LDAupper, LDAmean, LDAlower, Sign_Group, pvalue, fdr) %>% filter(!is.na(fdr))

#mp_plot_diff_res函数以分步呈现
#plot treeskeleton
p1 <- ggtree(taxa.tree,layout="radial",size = 0.3) +
  geom_point(
    data = td_filter(!isTip),
    fill="white",
    size=1,
    shape=21)
# display the high light of phylum clade.
p2 <- p1 +
  geom_hilight(mapping = aes(subset = nodeClass == "phylum", 
                             node = node, 
                             fill = label))
# display the relative abundance of features(OTU)
#按样本
p3 <- p2+
  ggnewscale::new_scale("fill") +
  geom_fruit(
    data = td_unnest(RareAbundanceBySample),
    geom = geom_star,
    mapping = aes(x = fct_reorder(Sample, Group, .fun=min),
                  size = RelRareAbundanceBySample,
                  fill = Group,
                  subset = RelRareAbundanceBySample > 0.05),
    starshape = 13,
    starstroke = 0.25,
    offset = 0.04,
    pwidth = 0.8,
    grid.params = list(linetype=2)) +
  scale_size_continuous(name="Relative Abundance (%)",range = c(1, 3)) +
  scale_fill_manual(values=c("#8b66b8","#c63596","#51c4c2","#2AB34A"))
p3
# display the tip labels of taxa tree
p4 <- p3 +
  geom_tiplab(size=2, offset=7.5)
p4
# display the LDA of significant OTU.
p5 <- p4 +
  ggnewscale::new_scale("fill") +
  geom_fruit(
    geom = geom_col,
    mapping = aes(x = LDAmean,
                  fill = Sign_Group),
    orientation = "y",
    offset = 0.4,
    pwidth = 0.5,
    axis.params = list(axis = "x",
                       title = "Log10(LDA)",
                       title.height = 0.01,
                       title.size = 2,
                       text.size = 1.8,
                       vjust = -1.5),
    grid.params = list(linetype = 2)
  )
p5
# display the significant (FDR) taxonomy after kruskal.test (default)
p6 <- p5 +
  ggnewscale::new_scale("size") +
  geom_point(
    data=td_filter(!is.na(fdr)),
    mapping = aes(size = -log10(pvalue),
                  fill = Sign_Group),
    shape = 21) +
  scale_size_continuous(range=c(1, 3)) +
  scale_fill_manual(values=c("#8b66b8","#51c4c2"))
p6 + theme(
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.3, "cm"),
  legend.spacing.y = unit(0.02, "cm"),
  legend.text = element_text(size = 7),
  legend.title = element_text(size = 9),
)

#修改p3和p4，其余不变
# display the relative abundance of features(OTU)
# 按分组
p3 <- p2+
  ggnewscale::new_scale("fill") +
  geom_fruit(
    data = td_unnest(RareAbundanceByGroup),
    geom = geom_star,
    mapping = aes(x=Group,
                  size = RelRareAbundanceByGroup,
                  fill = Group,
                  subset = RelRareAbundanceByGroup> 0.05),
    starshape = 13,
    starstroke = 0.25,
    offset = 0.04,
    pwidth = 0.1,
    grid.params = list(linetype=2)) +
  scale_size_continuous(name="Relative Abundance (%)",range = c(1, 3)) +
  scale_fill_manual(values=c("#8b66b8","#c63596","#51c4c2","#2AB34A"))
p3
# display the tip labels of taxa tree
p4 <- p3 +
  ggnewscale::new_scale("size") +
  geom_point(
    data=td_filter(!is.na(fdr)),
    mapping = aes(size = -log10(pvalue),
                  fill = Sign_Group),
    shape = 21) +
  scale_size_continuous(range=c(1, 3)) +
  scale_fill_manual(values=c("#8b66b8","#c63596","#51c4c2","#2AB34A"))
p4

p5 <- p4 +
  ggnewscale::new_scale("fill") +
  geom_fruit(
    geom = geom_col,
    mapping = aes(x = LDAmean,
                  fill = Sign_Group),
    orientation = "y",
    offset = 0.4,
    pwidth = 0.5,
    axis.params = list(axis = "x",
                       title = "Log10(LDA)",
                       title.height = 0.01,
                       title.size = 2,
                       text.size = 1.8,
                       vjust = -1.5),
    grid.params = list(linetype = 2)
  )
p5

p6 <- p5 +
  geom_tiplab(size=2, offset=2,mapping = aes(subset= LDAmean>=3))
p6
ggsave("C:/graduate/WH赤霉病/pdf picture/16S LefSe(2).pdf", p6,
       height=10,width =12,limitsize = FALSE )

######################Fungal LefSe########################3
setwd("C:/graduate/WH赤霉病/扩增子/WH Fungi")

rm(list=ls())
otumat <- read.csv('da.rar.csv', row.names = 1)
otumat<-as.data.frame(t(otumat))
otumat = data.matrix(otumat)

taxmat <- read.csv('LEfSe真菌ID.csv', row.names = 1)
taxmat = as.matrix(taxmat)

sampledata <- read.csv('真菌LEfSe分组信息.csv', row.names = 1)
sampledata = sample_data(data.frame(sampledata))


OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(taxmat)
OTU
TAX
sampledata
#加入phyloseq体系
physeq = phyloseq(OTU, TAX,sampledata)
physeq


set.seed(1024)#由于LDA效应大小是通过随机重采样计算的，因此应设置随机种子以实现结果的重复性
deres <- diff_analysis(obj = physeq, 
                       classgroup = "Group",#分组
                       firstcomfun = "kruskal_test",
                       padjust = "fdr",#p值校正方法
                       filtermod = "pvalue",#以pvalue列过滤
                       firstalpha = 0.05,
                       strictmod = TRUE,#是否进行一对一事后检验
                       secondcomfun = "wilcox_test",
                       secondalpha = 0.01,
                       mlfun = "lda",#线性判别分析，可选随机森林
                       ldascore=3#线性判别分数
)
deres


p<-ggdiffclade(obj=deres,
            layout="circular",#布局类型
            alpha=0.3, #树分支的背景透明度
            linewd=0.2, #树中连线粗细
            skpointsize=0.8, #树骨架中国点的大小
            taxlevel=2, #要展示的树分支水平2（门）及3以下（纲目科属种）
            cladetext=4,#文本大小
            setColors=F,#自定义颜色
            removeUnknown=T,#不删分支，但移除分类中有un_的物种注释
            reduce=T)+ # 移除分类不明的物种分支和其注释，为T时removeUnknown参数无效
  scale_fill_manual(values=c("#8b66b8","#51c4c2")) +
  guides(color = guide_legend(keywidth = 0.1, keyheight = 0.6,order = 3,ncol=1)) +
  theme(panel.background=element_rect(fill=NA),
        legend.position="right", 
        plot.margin=margin(0,0,0,0),
        legend.spacing.y=unit(0.1, "cm"), 
        legend.title=element_text(size=20),
        legend.text=element_text(size=16), 
        legend.box.spacing=unit(0.1,"cm"))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS LefSe(1).pdf", p,
       height=8,width =12,limitsize = FALSE )

# visualization of different results by ggdiffbox
diffbox <- ggdiffbox(obj=deres, box_notch=FALSE, 
                     colorlist=c("#8b66b8","#c63596","#51c4c2","#2AB34A"), 
                     l_xlabtext="Relative abundance")
diffbox
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS LefSe biomarker(1).pdf", diffbox,
       height=8,width =9,limitsize = FALSE )

effectsize <- ggeffectsize(obj=deres, lineheight=0.1,linewidth=0.3,pointsize=3) + 
  scale_color_manual(values=c("#8b66b8","#51c4c2")) 
effectsize
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS LefSe biomarker(2).pdf", effectsize,
       height=8,width =9,limitsize = FALSE )

ggdifftaxbar(obj=deres, xtextsize=1.5,figwidth = 10, ylabel = "Relative abundance",output="each_biomarker_barplot",coloslist=c("#8b66b8","#c63596","#51c4c2","#2AB34A"))

sss <- as.MPSE(physeq)
sss %<>% mp_rrarefy() 
sss %<>% mp_cal_abundance(.abundance = RareAbundance,force=T) %>% mp_cal_abundance(.abundance=RareAbundance,.group=Group,force = T)
sss %<>% mp_diff_analysis(
  .abundance = RelRareAbundanceBySample,
  .group = Group,
  p.adjust="fdr",
  filter.p="pvalue",
  first.test.alpha = 0.05,
  second.test.alpha = 0.01,
  ldascore = 3
)
taxa.tree <- sss %>% mp_extract_tree(type="taxatree")
taxa.tree %>% select(label, nodeClass, LDAupper, LDAmean, LDAlower, Sign_Group, pvalue, fdr) %>% filter(!is.na(fdr))

#mp_plot_diff_res函数以分步呈现
#plot treeskeleton
p1 <- ggtree(taxa.tree,layout="radial",size = 0.3) +
  geom_point(
    data = td_filter(!isTip),
    fill="white",
    size=1,
    shape=21)
# display the high light of phylum clade.
p2 <- p1 +
  geom_hilight(mapping = aes(subset = nodeClass == "phylum", 
                             node = node, 
                             fill = label))
# display the relative abundance of features(OTU)
#按样本
p3 <- p2+
  ggnewscale::new_scale("fill") +
  geom_fruit(
    data = td_unnest(RareAbundanceBySample),
    geom = geom_star,
    mapping = aes(x = fct_reorder(Sample, Group, .fun=min),
                  size = RelRareAbundanceBySample,
                  fill = Group,
                  subset = RelRareAbundanceBySample > 0.05),
    starshape = 13,
    starstroke = 0.25,
    offset = 0.04,
    pwidth = 0.8,
    grid.params = list(linetype=2)) +
  scale_size_continuous(name="Relative Abundance (%)",range = c(1, 3)) +
  scale_fill_manual(values=c("#8b66b8","#c63596","#51c4c2","#2AB34A"))
p3
# display the tip labels of taxa tree
p4 <- p3 +
  geom_tiplab(size=2, offset=7.5)
p4
# display the LDA of significant OTU.
p5 <- p4 +
  ggnewscale::new_scale("fill") +
  geom_fruit(
    geom = geom_col,
    mapping = aes(x = LDAmean,
                  fill = Sign_Group),
    orientation = "y",
    offset = 0.4,
    pwidth = 0.5,
    axis.params = list(axis = "x",
                       title = "Log10(LDA)",
                       title.height = 0.01,
                       title.size = 2,
                       text.size = 1.8,
                       vjust = -1.5),
    grid.params = list(linetype = 2)
  )
p5
# display the significant (FDR) taxonomy after kruskal.test (default)
p6 <- p5 +
  ggnewscale::new_scale("size") +
  geom_point(
    data=td_filter(!is.na(fdr)),
    mapping = aes(size = -log10(pvalue),
                  fill = Sign_Group),
    shape = 21) +
  scale_size_continuous(range=c(1, 3)) +
  scale_fill_manual(values=c("#8b66b8","#51c4c2"))
p6 + theme(
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.3, "cm"),
  legend.spacing.y = unit(0.02, "cm"),
  legend.text = element_text(size = 7),
  legend.title = element_text(size = 9),
)

#修改p3和p4，其余不变
# display the relative abundance of features(OTU)
# 按分组
p3 <- p2+
  ggnewscale::new_scale("fill") +
  geom_fruit(
    data = td_unnest(RareAbundanceByGroup),
    geom = geom_star,
    mapping = aes(x=Group,
                  size = RelRareAbundanceByGroup,
                  fill = Group,
                  subset = RelRareAbundanceByGroup> 0.05),
    starshape = 13,
    starstroke = 0.25,
    offset = 0.04,
    pwidth = 0.1,
    grid.params = list(linetype=2)) +
  scale_size_continuous(name="Relative Abundance (%)",range = c(1, 3)) +
  scale_fill_manual(values=c("#8b66b8","#c63596","#51c4c2","#2AB34A"))
p3
# display the tip labels of taxa tree
p4 <- p3 +
  ggnewscale::new_scale("size") +
  geom_point(
    data=td_filter(!is.na(fdr)),
    mapping = aes(size = -log10(pvalue),
                  fill = Sign_Group),
    shape = 21) +
  scale_size_continuous(range=c(1, 3)) +
  scale_fill_manual(values=c("#8b66b8","#c63596","#51c4c2","#2AB34A"))
p4

p5 <- p4 +
  ggnewscale::new_scale("fill") +
  geom_fruit(
    geom = geom_col,
    mapping = aes(x = LDAmean,
                  fill = Sign_Group),
    orientation = "y",
    offset = 0.4,
    pwidth = 0.5,
    axis.params = list(axis = "x",
                       title = "Log10(LDA)",
                       title.height = 0.01,
                       title.size = 2,
                       text.size = 1.8,
                       vjust = -1.5),
    grid.params = list(linetype = 2)
  )
p5

p6 <- p5 +
  geom_tiplab(size=2, offset=2,mapping = aes(subset= LDAmean>=3))
p6
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS LefSe(2).pdf", p6,
       height=10,width =12,limitsize = FALSE )
