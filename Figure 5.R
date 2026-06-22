library(ecodist)
library(pheatmap)
library(vegan)


setwd("C:/graduate/WH赤霉病/宏基因组/Metagene2-WH/Metagene2-WH/物种与功能")

#####################################################
rm(list=ls())
data2 <- read.csv("KO-tax.csv",header = T,row.names = 1)


#c查看其标准差
apply(data2,1,sd)

# 删掉标准差为0的行
data2 = data2[apply(data2, 1, function(x) sd(x)!=0),] 

# 删掉标准差为0的列
data2 = data2[,apply(data2, 2, function(x) sd(x)!=0)] 

mycolors <- colorRampPalette(c('#053061', '#68A8CF', 'white', '#EDB8B0', "red3"))(50)

tax_group <- read.delim('tax_group.txt', sep = '\t', row.names = 1)
gene_group <-read.delim("gene_group.txt",sep = '\t', row.names = 1)

ann_colors = list(Network=c("in network"="lightskyblue","not in network"="lightgoldenrod1"),
                  
                  Category =c("2, 3-butanediol"="green4","Siderophore"="orange2","Endochitinases"="purple2",
                              "Lysozyme"="dodgerblue2","Exochitinases"="yellow","Chitosanase"="lightpink","Chitin deacetylase"="lightcoral","Iturin A"="violet",
                              "Amphotericin B"="mediumpurple1","Syringomycin"="darkseagreen3","Surfactin"="seagreen1","Fengycin"="lightpink2","Endochitinases/Lysozyme"="paleturquoise1"
                  ))


#####聚类
drows<-distance(data2, "bray-curtis")


pdf("heatmap KO_tax.pdf", height =20, width = 25)
pheatmap(scale="row",as.matrix(data2),fontsize=10,fontsize_row=10,cex=1,cellwidth=10,cellheight =10,col=mycolors,annotation_row = gene_group, 
         annotation_col=tax_group,
         annotation_colors = ann_colors,
         border=FALSE,show_rownames=T,cluster_cols=F,cluster_rows=F,show_colnames = T)
dev.off()
