#################Cross-kingdom network####################
library(vegan)
library(psych)##
library(igraph)
library(beepr)

setwd("C:/graduate/WH赤霉病/宏基因组/kraken")
rm(list=ls())

da.rar<-read.csv("Bac-Fun.da.rar2.csv",head = T, row.names = 1)
env<-read.csv("env.csv",head = T, row.names = 1)
BF0<-read.csv("Bac-Fun.ID.csv", head = T, row.names = 1)
ID.tmp<-BF0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus", "species","OTU.ID")]

d1.raw<-da.rar[env$State=="Symptomatic",]
d2.raw<-da.rar[env$State=="Symptomless",]

fq <- 3; abu <- 8432; 

d1 <- d1.raw[,specnumber(t(d1.raw)) > fq & colSums(d1.raw)> abu]
d2 <- d2.raw[,specnumber(t(d2.raw)) > fq & colSums(d2.raw)> abu]

spman.d1 = corr.test(d1, use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
spman.d2 = corr.test(d2, use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)

######################################
# significant positive and negative correlations #
#####################################
BF0<-read.csv("Bac-Fun.ID.csv", head = T)
ID.tmp<-BF0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus","species","OTU.ID")]
r.cutoff = 0.6
p.cutoff = 0.05

spman.r0 <- spman.d1
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[lower.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[lower.tri(Cor)]], Cor=Cor[lower.tri(Cor)])

P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[lower.tri(P0)]], 
                 col=colnames(P0)[col(P0)[lower.tri(P0)]], p=P0[lower.tri(P0)])

df <- data.frame(Cor.df,  P.df, Compartment = "Rhizosphere", State = "Symptomatic")
da.tmp<-df.sig<- df[(df$Cor > r.cutoff | df$Cor < -r.cutoff) & df$p < p.cutoff,] 
da.tmp$color[da.tmp$Cor>r.cutoff]<-"darkorange1"
da.tmp$color[da.tmp$Cor< -r.cutoff]<-"green"
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp1<-data.frame(da.tmp, M1, M2)
write.csv(df.tmp1, "Diseaed Rhizosphere Bac-Fun Fusarium（正负）.csv")

spman.r0 <- spman.d2
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[lower.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[lower.tri(Cor)]], Cor=Cor[lower.tri(Cor)])

P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[lower.tri(P0)]], 
                 col=colnames(P0)[col(P0)[lower.tri(P0)]], p=P0[lower.tri(P0)])

df <- data.frame(Cor.df,  P.df, Compartment = "Rhizosphere", State = "Symptomless")
da.tmp<-df.sig<- df[(df$Cor > r.cutoff | df$Cor < -r.cutoff) & df$p < p.cutoff,] 
da.tmp$color[da.tmp$Cor>r.cutoff]<-"darkorange1"
da.tmp$color[da.tmp$Cor< -r.cutoff]<-"green"
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp2<-data.frame(da.tmp, M1, M2)
write.csv(df.tmp2, "Healthy Rhizosphere Bac-Fun Fusarium（正负）.csv")


####################
######bac-fun###
####################  

BF0<-read.csv("Bac-Fun network.ID（正负）.csv", head = T, row.names = 1)
ID.tmp<-BF0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus",	"OTU.ID","Group")]
ID.tmp$color <- "grey"
ID.tmp$color[ID.tmp$Group=="Fusarium"]<-"indianred1"
ID.tmp$color[ID.tmp$Group=="Bacteria"]<-"#33395B"
ID.tmp$color[ID.tmp$Group=="Fungi"]<-"#8E2D30"

da0<- df.tmp1
da <- da0[(da0$kingdom == "Bacteria"|da0$kingdom == "Fungi") & (da0$kingdom.1 == "Bacteria"|da0$kingdom.1 == "Fungi"),]
da[is.na(da)] = 0
da <-da[(da$kingdom!=0) & (da$kingdom.1!=0), ]

g <- graph.data.frame(da, directed=FALSE)
g.info = droplevels(ID.tmp[ID.tmp$OTU.ID %in% V(g)$name,])
g.info<-g.info[match(V(g)$name, g.info$OTU.ID),]
V(g)$color = as.character(g.info$color)
#V(g)$shape = as.character(g.info$shape)
write.graph(g ,"FR.bac-Fun_occor.gml",format="gml")

num.edges = length(E(g)) 
num.vertices = length(V(g))
connectance = edge_density(g,loops=FALSE)# 
degree<-igraph::degree(g)
degree1<-as.data.frame(degree)
average.degree = mean(igraph::degree(g))# 
average.path.length = average.path.length(g) 
diameter = diameter(g, directed = FALSE, unconnected = TRUE, weights = NULL)
edge.connectivity = edge_connectivity(g)
clustering.coefficient = transitivity(g) 
no.clusters = no.clusters(g)
centralization.betweenness = centralization.betweenness(g)$centralization 
centralization.degree = centralization.degree(g)$centralization
fun.fc<-cluster_fast_greedy(g)
Modularity<-modularity(fun.fc,membership(fun.fc))
No.modules<-nrow(data.frame(sizes(fun.fc)))
df.tmp1<-data.frame(network = "BF", "Rhizosphere", "Symptomatic",  num.edges, num.vertices, connectance, average.degree, average.path.length, diameter, edge.connectivity, clustering.coefficient,
                    no.clusters, centralization.betweenness,centralization.degree,  Modularity, No.modules)
write.csv(df.tmp1, "BF Diseased Rhizosphere Network（正负）.csv")

g.E <-data.frame(get.edgelist(g))
names(g.E)<- c("V1", "V2")

V1<-data.frame("v1"=g.E$V1)
V2<-data.frame("v2"=g.E$V2)

ID.tmpx<-ID.tmp[,c("OTU.ID", "kingdom")]
IDsub1<-ID.tmpx[ID.tmpx$OTU.ID %in% V1$v1, ]
IDsub2<-ID.tmpx[ID.tmpx$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]

M1$BF<-da$color

E(g)$color = as.character(M1$BF)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/宏基因组/kraken/cross-kingdom FR network.pdf",height=8,width =12)
plot(g, edge.width=1,  vertex.frame.color=NA,vertex.label=NA,edge.lty=1, edge.curved=T,vertex.size=3,main="Symptomatic Rhizosphere",cex.main=0.5) 
legend(title="Vertex",x=-1.8, y=0.2, 
       title.adj = 0,
       legend=c("Fusarium","Fungus(F)","Bacterium(B)"), 
       pch=c(21),
       pt.bg=c("indianred1","#8E2D30","#33395B"),
       col=c("indianred1","#8E2D30","#33395B"),
       pt.cex=1.5, cex=1, bty="n",box.lty=NA)
legend(title="Edge",x=-1.8, y=-0.2, title.col = "black",
       title.adj = 0,
       legend=c("Positive","Negative"),
       lty=1.5,col=c("darkorange1","green"),
       text.col = c("darkorange1","green"),
       cex=1, bty="n")
dev.off()

da0<- df.tmp2
da <- da0[(da0$kingdom == "Bacteria"|da0$kingdom == "Fungi") & (da0$kingdom.1 == "Bacteria"|da0$kingdom.1 == "Fungi"),]
da[is.na(da)] = 0
da <-da[(da$kingdom!=0) & (da$kingdom.1!=0), ]

g <- graph.data.frame(da, directed=FALSE)
g.info = droplevels(ID.tmp[ID.tmp$OTU.ID %in% V(g)$name,])
g.info<-g.info[match(V(g)$name, g.info$OTU.ID),]
V(g)$color = as.character(g.info$color)
#V(g)$shape = as.character(g.info$shape)
write.graph(g ,"HR.bac-Fun_occor.gml",format="gml")

num.edges = length(E(g)) 
num.vertices = length(V(g))
connectance = edge_density(g,loops=FALSE)# 
degree<-igraph::degree(g)
degree1<-as.data.frame(degree)
average.degree = mean(igraph::degree(g))# 
average.path.length = average.path.length(g) 
diameter = diameter(g, directed = FALSE, unconnected = TRUE, weights = NULL)
edge.connectivity = edge_connectivity(g)
clustering.coefficient = transitivity(g) 
no.clusters = no.clusters(g)
centralization.betweenness = centralization.betweenness(g)$centralization 
centralization.degree = centralization.degree(g)$centralization
fun.fc<-cluster_fast_greedy(g)
Modularity<-modularity(fun.fc,membership(fun.fc))
No.modules<-nrow(data.frame(sizes(fun.fc)))
df.tmp2<-data.frame(network = "BF", "Rhizosphere", "Symptomless",  num.edges, num.vertices, connectance, average.degree, average.path.length, diameter, edge.connectivity, clustering.coefficient,
                    no.clusters, centralization.betweenness,centralization.degree,  Modularity, No.modules)
write.csv(df.tmp2, "BF Healthy Rhizosphere Network（正负）.csv")

g.E <-data.frame(get.edgelist(g))
names(g.E)<- c("V1", "V2")

V1<-data.frame("v1"=g.E$V1)
V2<-data.frame("v2"=g.E$V2)

ID.tmpx<-ID.tmp[,c("OTU.ID", "kingdom")]
IDsub1<-ID.tmpx[ID.tmpx$OTU.ID %in% V1$v1, ]
IDsub2<-ID.tmpx[ID.tmpx$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]

M1$BF<-da$color

E(g)$color = as.character(M1$BF)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/宏基因组/kraken/cross-kingdom HR network.pdf",height=8,width =12)
plot(g, edge.width=1,  vertex.frame.color=NA,vertex.label=NA,edge.lty=1, edge.curved=T,vertex.size=3,main="Symptomless Rhizosphere",cex.main=0.5) 
legend(title="Vertex",x=-1.8, y=0.2, 
       title.adj = 0,
       legend=c("Fusarium","Fungus(F)","Bacterium(B)"), 
       pch=c(21),
       pt.bg=c("indianred1","#8E2D30","#33395B"),
       col=c("indianred1","#8E2D30","#33395B"),
       pt.cex=1.5, cex=1, bty="n",box.lty=NA)
legend(title="Edge",x=-1.8, y=-0.2, title.col = "black",
       title.adj = 0,
       legend=c("Positive","Negative"),
       lty=1.5,col=c("darkorange1","green"),
       text.col = c("darkorange1","green"),
       cex=1, bty="n")
dev.off()

################Network Pathway###############
library(ggplot2)

setwd("C:/graduate/WH赤霉病/宏基因组/kraken/1st_80")

rm(list=ls())
first80_ratio<-read.csv("1st80_ratio.csv",head = T)

p1 <- ggplot(first80_ratio, aes(x = reorder(Pathway, Abundance_ratio),y=Abundance_ratio)) +
  geom_col(fill="#0977B6",color="#0977B6",width = 0.3) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  labs(x = '', y = 'Abundance ratio(Symptomless/Symptomatic)')+
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x =element_text(color="black",size = 14,angle=0),
        axis.text.y =element_text(color="black",size = 14),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(face = "bold", colour = "black", size = 14),
        plot.title=element_text(face = "bold", colour = "black", size = 14))
p1
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/pathway picture/1st80_ratio.pdf", p1,
       height=4,width =12,limitsize = FALSE )

rm(list=ls())
first80_unique<-read.csv("1st80_unique.csv",head = T)

p2 <- ggplot(first80_unique, aes(x = reorder(Pathway, Abundance),y=Abundance)) +
  geom_col(aes(color=State,fill=State),width = 0.3) +
  scale_color_manual(values=c("#8b66b8","#51c4c2"))+
  scale_fill_manual(values=c("#8b66b8","#51c4c2"))+
  theme(legend.key=element_rect(colour=NA),panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  facet_grid(~State)+
  labs(x = '', y = 'Abundance')+
  theme(strip.text = element_text(size = 14,face="bold"),
        axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x =element_text(color="black",size = 14,angle=0),
        axis.text.y =element_text(color="black",size = 14),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(face = "bold", colour = "black", size = 14),
        plot.title=element_text(face = "bold", colour = "black", size = 14))
p2
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/pathway picture/1st80_unique.pdf", p2,
       height=12,width =15,limitsize = FALSE )

setwd("C:/graduate/WH赤霉病/宏基因组/kraken/2nd_80")
rm(list=ls())
second80_ratio<-read.csv("2nd80_ratio.csv",head = T)

p3 <- ggplot(second80_ratio, aes(x = reorder(Pathway, Abundance_ratio),y=Abundance_ratio)) +
  geom_col(fill="#0977B6",color="#0977B6",width = 0.3) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  labs(x = '', y = 'Abundance ratio(Symptomless/Symptomatic)')+
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x =element_text(color="black",size = 14,angle=0),
        axis.text.y =element_text(color="black",size = 14),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(face = "bold", colour = "black", size = 14),
        plot.title=element_text(face = "bold", colour = "black", size = 14))
p3
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/pathway picture/2nd80_ratio.pdf", p3,
       height=15,width =15,limitsize = FALSE )

rm(list=ls())
second80_unique<-read.csv("2nd80_unique.csv",head = T)

p4 <- ggplot(second80_unique, aes(x = reorder(Pathway, Abundance),y=Abundance)) +
  geom_col(aes(color=State,fill=State),width = 0.3) +
  scale_color_manual(values=c("#8b66b8","#51c4c2"))+
  scale_fill_manual(values=c("#8b66b8","#51c4c2"))+
  theme(legend.key=element_rect(colour=NA),panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  facet_grid(~State)+
  labs(x = '', y = 'Abundance')+
  theme(strip.text = element_text(size = 14,face="bold"),
        axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x =element_text(color="black",size = 14,angle=0),
        axis.text.y =element_text(color="black",size = 14),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(face = "bold", colour = "black", size = 14),
        plot.title=element_text(face = "bold", colour = "black", size = 14))
p4
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/pathway picture/2nd80_unique.pdf", p4,
       height=5,width =12,limitsize = FALSE )

setwd("C:/graduate/WH赤霉病/宏基因组/kraken/3rd_80")
rm(list=ls())
third80_ratio<-read.csv("3rd80_ratio.csv",head = T)

p5 <- ggplot(third80_ratio, aes(x = reorder(Pathway, Abundance_ratio),y=Abundance_ratio)) +
  geom_col(fill="#0977B6",color="#0977B6",width = 0.3) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  labs(x = '', y = 'Abundance ratio(Symptomless/Symptomatic)')+
  theme(axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x =element_text(color="black",size = 14,angle=0),
        axis.text.y =element_text(color="black",size = 14),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(face = "bold", colour = "black", size = 14),
        plot.title=element_text(face = "bold", colour = "black", size = 14))
p5
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/pathway picture/3rd80_ratio.pdf", p5,
       height=18,width =15,limitsize = FALSE )

rm(list=ls())
third80_unique<-read.csv("3rd80_unique.csv",head = T)

p6 <- ggplot(third80_unique, aes(x = reorder(Pathway, Abundance),y=Abundance)) +
  geom_col(aes(color=State,fill=State),width = 0.3) +
  scale_color_manual(values=c("#51c4c2"))+
  scale_fill_manual(values=c("#51c4c2"))+
  theme(legend.key=element_rect(colour=NA),panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  facet_grid(~State)+
  labs(x = '', y = 'Abundance')+
  theme(strip.text = element_text(size = 14,face="bold"),
        axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x =element_text(color="black",size = 14,angle=0),
        axis.text.y =element_text(color="black",size = 14),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(face = "bold", colour = "black", size = 14),
        plot.title=element_text(face = "bold", colour = "black", size = 14))
p6
ggsave("C:/graduate/WH赤霉病/宏基因组/kraken/pathway picture/3rd80_unique.pdf", p6,
       height=4,width =10,limitsize = FALSE )



##########################heatmap of KO and taxonomy#########################
library(ecodist)
library(pheatmap)
library(vegan)

setwd("C:/graduate/WH赤霉病/宏基因组/kraken/HR1st/KO_tax")

rm(list=ls())
data2 <- read.csv("KO_tax.csv",row.names = 1)


#c查看其标准差
apply(data2,1,sd)

# 删掉标准差为0的行
data2 = data2[apply(data2, 1, function(x) sd(x)!=0),] 

# 删掉标准差为0的列
data2 = data2[,apply(data2, 2, function(x) sd(x)!=0)] 

mycolors <- colorRampPalette(c('#053061', '#68A8CF', 'white', '#EDB8B0', "red3"))(50)


gene_group <-read.delim("KO_group.txt",sep = '\t', row.names = 1)

ann_colors = list(Category =c("UDP-N-acetylmuramoyl- pentapeptide biosynthesis I (meso-diaminopimelate containing)"="#FF0909",
                              "peptidoglycan biosynthesis I (meso-diaminopimelate containing)"="#31CD31",
                              "peptidoglycan maturation (meso-diaminopimelate containing)"="#FF77FF",
                              "peptidoglycan biosynthesis III (mycobacteria)"="#FFD912"
))


#####聚类
drows<-distance(data2, "bray-curtis")


pdf("heatmap KO_tax.pdf", height =6, width = 12)
pheatmap(scale="row",as.matrix(data2),fontsize=10,fontsize_row=10,cex=1,cellwidth=10,cellheight =10,col=mycolors,annotation_row = gene_group, 
         annotation_colors = ann_colors,
         border=FALSE,show_rownames=T,cluster_cols=F,cluster_rows=F,show_colnames = T)
dev.off()
