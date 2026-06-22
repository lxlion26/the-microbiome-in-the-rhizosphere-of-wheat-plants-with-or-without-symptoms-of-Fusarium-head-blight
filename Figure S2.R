#####################correlation###############
setwd("C:/graduate/WH赤霉病/宏基因组/从宏基因组中抽取16S")
library(vegan)
library(ggplot2)
library(ape)
library(splitstackshape)
library(reshape2)
library(dplyr)
library(linkET)

co<-read.csv("all.CO.csv",head=T,row.names = 1)
p1 <- qcorrplot(correlate(co, method = 'pearson'), type = 'lower', diag = FALSE) +  #环境变量矩阵计算 Spearman 相关系数
  geom_square() +  #绘制 Spearman 相关系数热图
  geom_mark(sep = '\n', size = 7, only_mark = T,sig_level = c(0.05, 0.01, 0.001),sig_thres = 0.05 ) +  #显示 Spearman 相关系数和显著性
  # geom_couple(aes(color = pd, size = rd), data = mantel, curvature = nice_curvature()) +  #环境和微生物的相关性展示为上述 Mantel 相关
  scale_fill_gradientn(colors = c('#053061', '#68A8CF', 'white', '#EDB8B0', '#C25759'), limits = c(-1, 1)) +  #根据 Spearman 相关指定热图颜色
  # scale_size_manual(values = c(0.5, 1, 2)) +  #根据 Mantel 相关指定线条粗细
  # scale_color_manual(values = c('#D95F02', '#1B9E77', '#E0E0E0')) +  #根据 Mantel 相关 p 值指定线条颜色
  guides(color = guide_legend(title = "Mantel's p", order = 1), #图例标题和排序
         size = guide_legend(title = "Mantel's r", order = 2), 
         fill = guide_colorbar(title = "Pearson's r", order = 3)) +
  theme(legend.key = element_blank())

p1
ggsave("C:/graduate/WH赤霉病/pdf picture/16S from metagenome correlation.pdf", p1,height=20,width =30,limitsize = FALSE )


############################Procrustes#################################
library(vegan)
library(ggplot2)

install.packages("C:/Users/shirl/Documents/R/win-library/4.1/SoDA_1.0-6.1.tar.gz",repos = NULL,type = "source")
library(SoDA)
library(picante)
#library(PCNM)
library(vegan)
library(ecodist)
library(MASS)
install.packages("C:/Users/shirl/Documents/R/win-library/4.1/flexmix_2.3-19.tar.gz",repos = NULL,type = "source")
install.packages("C:/Users/shirl/Documents/R/win-library/4.1/lmtest_0.9-40.tar.gz",repos = NULL,type = "source")
install.packages("C:/Users/shirl/Documents/R/win-library/4.1/betareg_3.1-4.tar.gz",repos = NULL,type = "source")
install.packages("C:/Users/shirl/Documents/R/win-library/4.1/hier.part_1.0-6.tar.gz",repos = NULL,type = "source")
library(hier.part)
library(car)
install.packages("C:/Users/shirl/Documents/R/win-library/4.1/MuMIn_1.43.17.tar.gz",repos = NULL,type = "source")
library(MuMIn)


library(vegan)
library(PMCMRplus)
library(MASS)
library(permute)
library(lattice)
library(car)
library(nlme)
library(mgcv)
library(labdsv)
library(sciplot)
library(spdep)
library(agricolae)

setwd("C:/graduate/WH赤霉病/宏基因组/从宏基因组中抽取16S")

###产生坐标数据
###############################################################################################################################################
#16S-OTU
rm(list = ls())

data0<-read.csv("amplicon.procust.csv",head = T, row.names = 1)
data0<-as.data.frame(t(data0))
data0<-data0[which(rowSums(data0) > 0),]
ID0<-read.csv("ID.onlyBac.扩增子.csv", head = T, row.names = 1)
env0<-read.csv("env.csv",head = T, row.names = 1)

###########   NMDS   #############
data.ra0<-data.frame(apply(data0,2,function(x) {x/sum(x)}))
data.ra<-as.data.frame(t(data.ra0))
data.ra.log0<-data.frame(apply(data.ra0,2,function(x) log(1+x)))#data stablize
#write.csv(data.ra.log,"data.ra.log.csv")
data.ra.log <- as.data.frame(t(data.ra.log0))

set.seed(999)

mds.run <- metaMDS(data.ra.log,trymax=100)
mds<-data.frame(mds.run$points)
mds<-tibble::rownames_to_column(mds,var = "Sample")
mds.env<-merge(mds,env0,by="Sample")
write.csv(mds.env,"amplicon.nmds.csv")##整理成txt格式

#############################################################################################################################################
#16S-KO
rm(list = ls())

data0<-read.csv("16S.procust.csv",head = T, row.names = 1)
data0<-as.data.frame(t(data0))
data0<-data0[which(rowSums(data0) > 0),]
ID0<-read.csv("ID.onlyBac.csv", head = T, row.names = 1)
env0<-read.csv("env.csv",head = T, row.names = 1)

###########   NMDS   #############
data.ra0<-data.frame(apply(data0,2,function(x) {x/sum(x)}))
data.ra<-as.data.frame(t(data.ra0))
data.ra.log0<-data.frame(apply(data.ra0,2,function(x) log(1+x)))#data stablize
#write.csv(data.ra.log,"data.ra.log.csv")
data.ra.log <- as.data.frame(t(data.ra.log0))

set.seed(999)

mds.run <- metaMDS(data.ra.log,trymax=100)
mds<-data.frame(mds.run$points)
mds<-tibble::rownames_to_column(mds,var = "Sample")
mds.env<-merge(mds,env0,by="Sample")
write.csv(mds.env,"16S.nmds.csv")##整理成txt格式

##################################################################正式开始分析#######################################################

##读入坐标数据
rm(list = ls())
A <- read.csv("amplicon.nmds.csv",head = T, row.names = 1)
R <- read.csv("16S.nmds.csv",head = T, row.names = 1)
A<-A[,c(1,2)]
R<-R[,c(1,2)]
#执行 Procrustes 分析，详情 ?procrustes
#以对称分析为例（symmetric = TRUE）
proc <- procrustes(Y = A, X = R, symmetric = TRUE)
summary(proc)

#旋转图
plot(proc, kind = 1, type = 'text')
#一些重要的结果提取
names(proc)

head(proc$Yrot)  #Procrustes 分析后 Y 的坐标
head(proc$X)  #Procrustes 分析后 X 的坐标
proc$ss  #偏差平方和 M2 统计量
proc$rotation  #通过该值可获得旋转轴的坐标位置
#残差图
plot(proc, kind = 2)
residuals(proc)  #残差值
#PROTEST 检验，详情 ?protest
#以 999 次置换为例
#注：protest() 中执行的是对称 Procrustes 分析，X 和 Y 的分配调换不影响 M2 统计量的计算
set.seed(123)
prot <- protest(Y = A, X = R, permutations = how(nperm = 999))
prot

#重要统计量的提取
names(prot)
prot$signif  #p 值
prot$ss  #偏差平方和 M2 统计量

#提取 Procrustes 分析的坐标
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)

#添加分组信息
group <- read.csv("group.csv",head = T)

Y$Samples <- rownames(Y)
Y <- merge(Y, group, by = 'Samples')
write.csv(Y,"Y.csv")##整理成txt格式,加入type信息
#Y <- as.data.frame(Y)

#####寻找边界
#df_points <- read.csv("Y.16S.KO.csv",head=T,row.names=1)
#library(plyr)
#df <-df_points[,c("X1", "X2", "Type")]
#find_hull <- function(df_points) df_points[chull(df_points$X1, df_points$X2), ]
#hulls1 <- ddply(df, "Type", find_hull )###寻找边界
#write.csv(hulls1,"hulls1.16S.KO.csv")


#df2 <-df_points[,c("x", "y", "Type")]
#find_hull2 <- function(df_points) df_points[chull(df_points$x, df_points$y), ]
#hulls2 <- ddply(df2, "Type", find_hull2 )###寻找边界
#write.csv(hulls2,"hulls2.16S.KO.csv")



####作图尝试
#hulls <- read.csv("hulls.16S.KO.csv",head=T,row.names=1)
Y <- read.csv("Y.csv",head=T,row.names=1)
Y2 <- read.csv("Y2.csv",head=T,row.names=1)
##边界
df_points <- read.csv("Y2.csv",head=T,row.names=1)
library(plyr)
df <-df_points[,c("X", "Y", "State")]
find_hull <- function(df_points) df_points[chull(df_points$X, df_points$Y), ]
hulls <- ddply(df, "State", find_hull )###寻找边界



#####作图

p1<-ggplot(data=Y2)+ 
  geom_point(aes(x=X, y=Y, color = State,shape=method), size = 3.5,alpha=0.7) +
  
  scale_color_manual(values = c("#8b66b8","#51c4c2"), 
                     limits = c('Symptomatic','Symptomless'))+
  scale_shape_manual(values=c(19, 17))+
  
  geom_polygon(data=hulls,aes(x=X,y=Y,fill = State), alpha=.1)+
  #geom_polygon(data=hulls2, alpha=.1,aes(fill = Type))+
  scale_fill_manual(values=c("#8b66b8","#51c4c2"))+
  
  
  geom_segment(data=Y, aes(x = X1, y = X2, xend = x, yend = y,color = State), 
               arrow = arrow(length = unit(0.3, 'cm')),size = 0.6) +
  
  theme_bw()+
  guides(color=guide_legend(title= "State"))+
  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent')) +
  
  theme(legend.position = 'right',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=20, face="bold"),
        legend.text = element_text(colour="black", size=16),
        axis.text=element_text(size=20,face="bold",colour="black"),
        axis.title=element_text(size=20,face="bold"),
        plot.title=element_text(face = "bold", colour = "black", size = 20))+
  labs(x = 'NMDS1', y = 'NMDS2', color = '') +
  geom_vline(xintercept = 0, color = 'gray34', linetype = 2, linewidth = 1) +
  geom_hline(yintercept = 0, color = 'gray34', linetype = 2, linewidth = 1) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], linewidth = 1.3) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], linewidth = 1.3) +
  
  ggtitle(expression(M ^2~'= 0.721***'~~~italic(P)~'= 0.001'))
p1
ggsave("C:/graduate/WH赤霉病/pdf picture/procrust 16S-amplicon.pdf", p1,height=8,width =10,limitsize = FALSE )
