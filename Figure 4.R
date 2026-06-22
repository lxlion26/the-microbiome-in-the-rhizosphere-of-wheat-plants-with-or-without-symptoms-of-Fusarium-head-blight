#################Cross-kingdom network####################
library(vegan)
library(psych)##
library(igraph)
library(beepr)

setwd("C:/graduate/WH赤霉病/扩增子")
rm(list=ls())

da.rar<-read.csv("Bac-Fun.da.rar.csv",head = T, row.names = 1)
env<-read.csv("Bac-Fun.env.csv",head = T, row.names = 1)
BF0<-read.csv("Bac-Fun.ID.csv", head = T, row.names = 1)
ID.tmp<-BF0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus", "OTU.ID")]


env.R<-env[env$Compartment=="Rhizosphere",]
env.Z<-env[env$Compartment=="Root",]
d1.R.raw<-da.rar[env$Compartment=="Rhizosphere",]
d2.Z.raw<-da.rar[env$Compartment=="Root",]

fq <- 4; abu <- 0; 

d1 <- d1.R.raw [env.R$State=="Symptomatic", ]
d1 <- d1 [,specnumber(t(d1)) > fq & colSums(d1)> abu]


d2 <- d1.R.raw [env.R$State=="Symptomless", ]
d2 <- d2 [,specnumber(t(d2)) > fq & colSums(d2)> abu]

d3 <- d2.Z.raw [env.Z$State=="Symptomatic", ]
d3 <- d3 [,specnumber(t(d3)) > fq & colSums(d3)> abu]


d4 <- d2.Z.raw [env.Z$State=="Symptomless", ]
d4 <- d4 [,specnumber(t(d4)) > fq & colSums(d4)> abu]


spman.d1 = corr.test(d1, use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
spman.d2 = corr.test(d2, use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
spman.d3 = corr.test(d3, use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)
spman.d4 = corr.test(d4, use="pairwise",method="spearman",adjust="fdr", alpha=.05, ci=FALSE)

######################################
# significant positive and negative correlations #
#####################################
BF0<-read.csv("Bac-Fun.ID.csv", head = T, row.names = 1)
ID.tmp<-BF0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus","OTU.ID")]
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
#write.csv(df.tmp1, "Diseaed Rhizosphere Bac-Fun Fusarium（正负）.csv")

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
#write.csv(df.tmp2, "Healthy Rhizosphere Bac-Fun Fusarium（正负）.csv")

spman.r0 <- spman.d3
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[lower.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[lower.tri(Cor)]], Cor=Cor[lower.tri(Cor)])
P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[lower.tri(P0)]], 
                 col=colnames(P0)[col(P0)[lower.tri(P0)]], p=P0[lower.tri(P0)])

df <- data.frame(Cor.df,  P.df, Compartment = "Root", State = "Symptomatic")
da.tmp<-df.sig<- df[ (df$Cor > r.cutoff | df$Cor < -r.cutoff) & df$p < p.cutoff,] 
da.tmp$color[da.tmp$Cor>r.cutoff]<-"darkorange1"
da.tmp$color[da.tmp$Cor< -r.cutoff]<-"green"
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp3<-data.frame(da.tmp, M1, M2)
#write.csv(df.tmp3, "Diseased Root Bac-Fun Fusarium（正负）.csv")

spman.r0 <- spman.d4
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[lower.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[lower.tri(Cor)]], Cor=Cor[lower.tri(Cor)])

P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[lower.tri(P0)]], 
                 col=colnames(P0)[col(P0)[lower.tri(P0)]], p=P0[lower.tri(P0)])
df <- data.frame(Cor.df,  P.df, Compartment = "Root", State = "Symptomless")
da.tmp<-df.sig<- df[ (df$Cor > r.cutoff | df$Cor < -r.cutoff) & df$p < p.cutoff,] 
da.tmp$color[da.tmp$Cor>r.cutoff]<-"darkorange1"
da.tmp$color[da.tmp$Cor< -r.cutoff]<-"green"
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp4<-data.frame(da.tmp, M1, M2) 
#write.csv(df.tmp4, "Healthy Root Bac-Fun Fusarium（正负）.csv")

####################
######bac-fun###
####################   
BF0<-read.csv("Diseased Rhizosphere Bac-Fun.ID（正负）.csv", head = T, row.names = 1)
ID.tmp<-BF0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus",	"OTU.ID","Group")]
ID.tmp$shape[ID.tmp$kingdom=="Bacteria"]<-"circle"
ID.tmp$shape[ID.tmp$kingdom=="Fungi"]<-"square"
ID.tmp$color <- "grey"
ID.tmp$color[ID.tmp$genus=="Fusarium"]<-"indianred1"
ID.tmp$color[ID.tmp$Group=="Plant pathogens"]<-"#8b66b8"
ID.tmp$color[ID.tmp$Group=="Potential beneficial taxa"]<-"#51c4c2"

da0<- df.tmp1
da <- da0[(da0$kingdom == "Bacteria"|da0$kingdom == "Fungi") & (da0$kingdom.1 == "Bacteria"|da0$kingdom.1 == "Fungi"),]
da[is.na(da)] = 0
da <-da[(da$kingdom!=0) & (da$kingdom.1!=0), ]


g <- graph.data.frame(da, directed=FALSE)
g.info = droplevels(ID.tmp[ID.tmp$OTU.ID %in% V(g)$name,])
g.info<-g.info[match(V(g)$name, g.info$OTU.ID),]
V(g)$color = as.character(g.info$color)
V(g)$shape = as.character(g.info$shape)

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
pdf("C:/graduate/WH赤霉病/pdf picture/cross-kingdom FR network.pdf",height=8,width =12)
plot(g, edge.width=1,  vertex.frame.color=NA,vertex.label=NA,edge.lty=1, edge.curved=T,vertex.size=3,main="Symptomatic Rhizosphere",cex.main=0.5) 
legend(title="Vertex color",x=-1.8, y=0.5, 
       title.adj = 0,
       legend=c("Fusarium","Plant pathogens","Potential beneficial taxa"),
       pt.bg=c("indianred1","#8b66b8","#51c4c2"),
       col=c("indianred1","#8b66b8","#51c4c2"), 
       pch=24,pt.cex=1.5, cex=1, bty="n",box.lty=NA)
legend(title="Vertex shape",x=-1.8, y=0.1, 
       title.adj = 0,
       legend=c("Fungus(F)","Bacterium(B)"), 
       pch=c(22,21),
       col=c("grey"), 
       pt.bg=c("grey"),
       pt.cex=1.5, cex=1, bty="n",box.lty=NA)
legend(title="Edge color",x=-1.8, y=-0.2, title.col = "black",
       title.adj = 0,
       legend=c("Positive","Negative"),
       lty=1.5,col=c("darkorange1","green"),
       text.col = c("darkorange1","green"),
       cex=1, bty="n")
dev.off()

BF0<-read.csv("Healthy Rhizosphere Bac-Fun.ID（正负）.csv", head = T, row.names = 1)
ID.tmp<-BF0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus",	"OTU.ID","Group")]
ID.tmp$shape[ID.tmp$kingdom=="Bacteria"]<-"circle"
ID.tmp$shape[ID.tmp$kingdom=="Fungi"]<-"square"
ID.tmp$color <- "grey"
ID.tmp$color[ID.tmp$genus=="Fusarium"]<-"indianred1"
ID.tmp$color[ID.tmp$Group=="Plant pathogens"]<-"#8b66b8"
ID.tmp$color[ID.tmp$Group=="Potential beneficial taxa"]<-"#51c4c2"

da0<- df.tmp2
da <- da0[(da0$kingdom == "Bacteria"|da0$kingdom == "Fungi") & (da0$kingdom.1 == "Bacteria"|da0$kingdom.1 == "Fungi"),]
da[is.na(da)] = 0
da <-da[(da$kingdom!=0) & (da$kingdom.1!=0), ]


g <- graph.data.frame(da, directed=FALSE)
g.info = droplevels(ID.tmp[ID.tmp$OTU.ID %in% V(g)$name,])
g.info<-g.info[match(V(g)$name, g.info$OTU.ID),]
V(g)$color = as.character(g.info$color)
V(g)$shape = as.character(g.info$shape)

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
pdf("C:/graduate/WH赤霉病/pdf picture/cross-kingdom HR network.pdf",height=8,width =12)
plot(g, edge.width=1,  vertex.frame.color=NA,vertex.label=NA,edge.lty=1, edge.curved=T,vertex.size=3,main="Symptomless Rhizosphere",cex.main=0.5) 
legend(title="Vertex color",x=-1.8, y=0.5, 
       title.adj = 0,
       legend=c("Fusarium","Plant pathogens","Potential beneficial taxa"),
       pt.bg=c("indianred1","#8b66b8","#51c4c2"),
       col=c("indianred1","#8b66b8","#51c4c2"), 
       pch=24,pt.cex=1.5, cex=1, bty="n",box.lty=NA)
legend(title="Vertex shape",x=-1.8, y=0.1, 
       title.adj = 0,
       legend=c("Fungus(F)","Bacterium(B)"), 
       pch=c(22,21),
       col=c("grey"), 
       pt.bg=c("grey"),
       pt.cex=1.5, cex=1, bty="n",box.lty=NA)
legend(title="Edge color",x=-1.8, y=-0.2, title.col = "black",
       title.adj = 0,
       legend=c("Positive","Negative"),
       lty=1.5,col=c("darkorange1","green"),
       text.col = c("darkorange1","green"),
       cex=1, bty="n")
dev.off()

BF0<-read.csv("Diseased Root Bac-Fun.ID（正负）.csv", head = T, row.names = 1)
ID.tmp<-BF0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus",	"OTU.ID","Group")]
ID.tmp$shape[ID.tmp$kingdom=="Bacteria"]<-"circle"
ID.tmp$shape[ID.tmp$kingdom=="Fungi"]<-"square"
ID.tmp$color <- "grey"
ID.tmp$color[ID.tmp$genus=="Fusarium"]<-"indianred1"
ID.tmp$color[ID.tmp$Group=="Plant pathogens"]<-"#8b66b8"
ID.tmp$color[ID.tmp$Group=="Potential beneficial taxa"]<-"#51c4c2"

da0<- df.tmp3
da <- da0[(da0$kingdom == "Bacteria"|da0$kingdom == "Fungi") & (da0$kingdom.1 == "Bacteria"|da0$kingdom.1 == "Fungi"),]
da[is.na(da)] = 0
da <-da[(da$kingdom!=0) & (da$kingdom.1!=0), ]


g <- graph.data.frame(da, directed=FALSE)
g.info = droplevels(ID.tmp[ID.tmp$OTU.ID %in% V(g)$name,])
g.info<-g.info[match(V(g)$name, g.info$OTU.ID),]
V(g)$color = as.character(g.info$color)
V(g)$shape = as.character(g.info$shape)

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
df.tmp3<-data.frame(network = "BF", "Root", "Symptomatic",  num.edges, num.vertices, connectance, average.degree, average.path.length, diameter, edge.connectivity, clustering.coefficient,
                    no.clusters, centralization.betweenness,centralization.degree,  Modularity, No.modules)
write.csv(df.tmp3, "BF Diseased Root Network（正负）.csv")

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
pdf("C:/graduate/WH赤霉病/pdf picture/cross-kingdom FZ network.pdf",height=8,width =12)
plot(g, edge.width=1,  vertex.frame.color=NA,vertex.label=NA,edge.lty=1, edge.curved=T,vertex.size=3,main="Symptomatic Root",cex.main=0.5) 
legend(title="Vertex color",x=-1.8, y=0.5, 
       title.adj = 0,
       legend=c("Fusarium","Plant pathogens","Potential beneficial taxa"),
       pt.bg=c("indianred1","#8b66b8","#51c4c2"),
       col=c("indianred1","#8b66b8","#51c4c2"), 
       pch=24,pt.cex=1.5, cex=1, bty="n",box.lty=NA)
legend(title="Vertex shape",x=-1.8, y=0.1, 
       title.adj = 0,
       legend=c("Fungus(F)","Bacterium(B)"), 
       pch=c(22,21),
       col=c("grey"), 
       pt.bg=c("grey"),
       pt.cex=1.5, cex=1, bty="n",box.lty=NA)
legend(title="Edge color",x=-1.8, y=-0.2, title.col = "black",
       title.adj = 0,
       legend=c("Positive","Negative"),
       lty=1.5,col=c("darkorange1","green"),
       text.col = c("darkorange1","green"),
       cex=1, bty="n")
dev.off()

BF0<-read.csv("Healthy Root Bac-Fun.ID（正负）.csv", head = T, row.names = 1)
ID.tmp<-BF0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus",	"OTU.ID","Group")]
ID.tmp$shape[ID.tmp$kingdom=="Bacteria"]<-"circle"
ID.tmp$shape[ID.tmp$kingdom=="Fungi"]<-"square"
ID.tmp$color <- "grey"
ID.tmp$color[ID.tmp$genus=="Fusarium"]<-"indianred1"
ID.tmp$color[ID.tmp$Group=="Plant pathogens"]<-"#8b66b8"
ID.tmp$color[ID.tmp$Group=="Potential beneficial taxa"]<-"#51c4c2"

da0<- df.tmp4
da <- da0[(da0$kingdom == "Bacteria"|da0$kingdom == "Fungi") & (da0$kingdom.1 == "Bacteria"|da0$kingdom.1 == "Fungi"),]
da[is.na(da)] = 0
da <-da[(da$kingdom!=0) & (da$kingdom.1!=0), ]


g <- graph.data.frame(da, directed=FALSE)
g.info = droplevels(ID.tmp[ID.tmp$OTU.ID %in% V(g)$name,])
g.info<-g.info[match(V(g)$name, g.info$OTU.ID),]
V(g)$color = as.character(g.info$color)
V(g)$shape = as.character(g.info$shape)

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
df.tmp4<-data.frame(network = "BF", "Root", "Symptomless",  num.edges, num.vertices, connectance, average.degree, average.path.length, diameter, edge.connectivity, clustering.coefficient,
                    no.clusters, centralization.betweenness,centralization.degree,  Modularity, No.modules)
write.csv(df.tmp4, "BF Healthy Root Network（正负）.csv")


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
pdf("C:/graduate/WH赤霉病/pdf picture/cross-kingdom HZ network.pdf",height=8,width =12)
plot(g, edge.width=1,  vertex.frame.color=NA,vertex.label=NA,edge.lty=1, edge.curved=T,vertex.size=3,main="Symptomless Root",cex.main=0.5) 
legend(title="Vertex color",x=-1.8, y=0.5, 
       title.adj = 0,
       legend=c("Fusarium","Plant pathogens","Potential beneficial taxa"),
       pt.bg=c("indianred1","#8b66b8","#51c4c2"),
       col=c("indianred1","#8b66b8","#51c4c2"), 
       pch=24,pt.cex=1.5, cex=1, bty="n",box.lty=NA)
legend(title="Vertex shape",x=-1.8, y=0.1, 
       title.adj = 0,
       legend=c("Fungus(F)","Bacterium(B)"), 
       pch=c(22,21),
       col=c("grey"), 
       pt.bg=c("grey"),
       pt.cex=1.5, cex=1, bty="n",box.lty=NA)
legend(title="Edge color",x=-1.8, y=-0.2, title.col = "black",
       title.adj = 0,
       legend=c("Positive","Negative"),
       lty=1.5,col=c("darkorange1","green"),
       text.col = c("darkorange1","green"),
       cex=1, bty="n")
dev.off()
########################Cross-kingdom network no.vertices and edges################
library(vegan)
library(ggplot2)
library(colorRamps)
library(ape)
library(splitstackshape)
library(reshape2)

setwd("C:/graduate/WH赤霉病/扩增子")
rm(list=ls())

NP <- read.csv("Fusarium network正负.csv")
col11<-c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#97D1A0", "#43A743", "grey","#ff00ff","#00ff00", "deepskyblue", "gold", "red", "navy", "darkgreen","maroon3", "black", "bisque", "grey")
p1<-ggplot(NP, aes(x = State, y = Vertex, fill=Group)) +
  geom_bar(stat='identity', position = "stack", width=0.5)+
  scale_fill_manual(values= col11)+
  theme_bw()+
  facet_grid(~Compartment)+
  labs(x="",y="Number of vertices")+
  guides(fill=guide_legend(title= "Group"))+
  #scale_y_continuous(limit = c(0, 1))+
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 20,face="bold"),
        legend.title = element_text(colour="black", size=20, face="bold"),
        legend.text = element_text(colour="black", size=16),
        axis.text.y=element_text(colour="black",size=16),
        axis.text.x=element_text(colour="black",size=16,angle = -20),
        axis.title=element_text(colour="black",size=20,face="bold"))
p1
ggsave("C:/graduate/WH赤霉病/pdf picture/cross-kingdom network num.vertices.pdf", p1,
       height=8,width =9,limitsize = FALSE )

p2<-ggplot(NP, aes(x = State, y = Edge, fill=Group)) +
  geom_bar(stat='identity', position = "stack", width=0.5)+
  scale_fill_manual(values= col11)+
  theme_bw()+
  facet_grid(~Compartment)+
  labs(x="",y="Number of edges")+
  guides(fill=guide_legend(title= "Group"))+
  #scale_y_continuous(limit = c(0, 1))+
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 20,face="bold"),
        legend.title = element_text(colour="black", size=20, face="bold"),
        legend.text = element_text(colour="black", size=16),
        axis.text.y=element_text(colour="black",size=16),
        axis.text.x=element_text(colour="black",size=16,angle = -20),
        axis.title=element_text(colour="black",size=20,face="bold"))
p2
ggsave("C:/graduate/WH赤霉病/pdf picture/cross-kingdom network num.edges.pdf", p2,
       height=8,width =9,limitsize = FALSE )

###########################ZIPI################################
library(vegan)##decostand
library(WGCNA)##corAndPvalue
library(multtest)##mt.rawp2adjp
library(igraph)##graph.adjacency    delete.vertices
library(brainGraph)##Zi-Pi
library(ggplot2)

setwd("C:/graduate/WH赤霉病/扩增子")
rm(list=ls())
bac.coco<-read.csv("FR.bac-fun.relabu.csv",head=T,row.names = 1)

matrix2igraph2<-function(matr,r.threshold,p.threshold){
  
  occor<-corAndPvalue(matr,method = c( "spearman")) 
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="fdr") 
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(matr)[2])
  occor.r<-occor$cor
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold]<-0 
  diag(occor.r)<-0  
  igraph<-graph.adjacency(occor.r,mode="undirected", weighted=TRUE, diag=FALSE)# NOTE:????????weighted=NULL,???Ǵ?ʱҪע???˺???ֻ??ʶ???໥???þ?????????????????Ӧ??ǰ??ȷ????????ȷ
  igraph<- igraph::simplify(igraph)
  igraph
}

bac<-matrix2igraph2(bac.coco,0.60,0.05)

bad.vs<-V(bac)[degree(bac) == 0] 
igraph <- delete.vertices(bac, bad.vs) 
igraph

##modularity 
E(igraph)$weight<-abs(E(igraph)$weight)
bac.fc<-cluster_fast_greedy(igraph)
modularity(bac.fc)##>0.4 indicate modular structures
modularity(bac.fc,membership(bac.fc))
membership(bac.fc)
sizes(bac.fc)


###Calculate vertex within-module degree z-score:???????????ġ?��?ӵ?
bac.tax_0.01<-read.csv("FR.Bac-Fun.IDotutab.rar.csv",head=T,row.names=1)
bac.g.tax_0.01<-bac.tax_0.01[V(igraph)$name,]
bac.g.tax_0.01$phylum<-as.factor(bac.g.tax_0.01$phylum)
levels(bac.g.tax_0.01$phylum)
bac.g.tax<-droplevels(bac.g.tax_0.01$phylum)
bac.g.tax
levels(bac.g.tax)<-c("black","blue","dodgerblue","dodgerblue4","gray26",'red','orange')????????
  levels(bac.g.tax)

str(bac.g.tax_0.01)

library(brainGraph)

bac.comps <- membership(bac.fc)
bac_Pi<-part_coeff(igraph, membership(bac.fc))
bac.Zi<-within_module_deg_z_score(igraph, membership(bac.fc))
write.csv(bac_Pi,"FR.bac-fun_Pi.csv")
write.csv(bac.Zi,"FR.bac-fun_Zi.csv")

FR.BF.ZIPI<-read.csv("FR.Bac-Fun.ZIPI.csv",head=T,row.names = 1)
FR.BF.ZIPI$kingdom<-as.factor(FR.BF.ZIPI$kingdom)
levels(FR.BF.ZIPI$kingdom)
FR.BF.ZIPI<-droplevels(FR.BF.ZIPI$kingdom)
#levels(FR.BF.ZIPI)<-c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#43A743","#97D1A0","grey","orange", "red","cyan")
levels(FR.BF.ZIPI)

FR.BF.ZIPI<-read.csv("FR.Bac-Fun.ZIPI.csv",head=T,row.names = 1)

p1<-ggplot(FR.BF.ZIPI, aes(x = Pi, y = Zi, color=kingdom,shape=type)) +
  geom_point(aes(colour=kingdom,shape=type),size=5,alpha=0.5)+ 
  scale_color_manual(values = c("#ABDAEC","#43A743"))+
  scale_shape_manual(values = c(18,17,1))+
  labs(title="Cross-kingdom Network",subtitle="Symptomatic Rhizosphere",x="Among-module connectivity (Pi)",y = "Within-module connectivity (Zi)")+
  scale_y_continuous(limits = c(-3,5))+
  scale_x_continuous(limits = c(0,1))+
  guides(fill=guide_legend(title= "kingdom"),shape=guide_legend(title= "type"))+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 20,face="bold"),
        legend.title = element_text(colour="black", size=20, face="bold"),
        legend.text = element_text(colour="black", size=16),
        axis.text.y=element_text(colour="black",size=16),
        axis.text.x=element_text(colour="black",size=16,angle = 0),
        axis.title=element_text(colour="black",size=20,face="bold"),
        plot.title=element_text(face = "bold", colour = "black", size = 20),
        plot.subtitle=element_text(colour = "black", size = 16))+
  geom_vline(xintercept = 0.62, color = 'gray34', linetype = 1, linewidth = 1.5) +
  geom_hline(yintercept = 2.5, color = 'gray34', linetype = 1, linewidth = 1.5) 
p1
ggsave("C:/graduate/WH赤霉病/pdf picture/cross-kingdom network FR ZIPI.pdf", p1,
       height=8,width =9,limitsize = FALSE )

rm(list=ls())
bac.coco<-read.csv("HR.bac-fun.relabu.csv",head=T,row.names = 1)
matrix2igraph2<-function(matr,r.threshold,p.threshold){
  
  occor<-corAndPvalue(matr,method = c( "spearman")) 
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="fdr") 
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(matr)[2])
  occor.r<-occor$cor
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold]<-0 
  diag(occor.r)<-0  
  igraph<-graph.adjacency(occor.r,mode="undirected", weighted=TRUE, diag=FALSE)
  igraph<- igraph::simplify(igraph)
  igraph
}

bac<-matrix2igraph2(bac.coco,0.60,0.05) 


bad.vs<-V(bac)[degree(bac) == 0] 
igraph <- delete.vertices(bac, bad.vs) 
igraph

##modularity 
E(igraph)$weight<-abs(E(igraph)$weight)
bac.fc<-cluster_fast_greedy(igraph)
modularity(bac.fc)##>0.4 indicate modular structures
modularity(bac.fc,membership(bac.fc))
membership(bac.fc)
sizes(bac.fc)


###Calculate vertex within-module degree z-score:???????????ġ?��?ӵ?
bac.tax_0.01<-read.csv("HR.Bac-Fun.IDotutab.rar.csv",head=T,row.names=1)
bac.g.tax_0.01<-bac.tax_0.01[V(igraph)$name,]
bac.g.tax_0.01$phylum<-as.factor(bac.g.tax_0.01$phylum)
levels(bac.g.tax_0.01$phylum)
bac.g.tax<-droplevels(bac.g.tax_0.01$phylum)
bac.g.tax
levels(bac.g.tax)<-c("black","blue","dodgerblue","dodgerblue4","gray26",'red','orange')
  levels(bac.g.tax)

str(bac.g.tax_0.01)

library(brainGraph)

bac.comps <- membership(bac.fc)
bac_Pi<-part_coeff(igraph, membership(bac.fc))
bac.Zi<-within_module_deg_z_score(igraph, membership(bac.fc))
write.csv(bac_Pi,"HR.bac-fun_Pi.csv")
write.csv(bac.Zi,"HR.bac-fun_Zi.csv")

HR.BF.ZIPI<-read.csv("HR.Bac-Fun.ZIPI.csv",head=T,row.names = 1)
HR.BF.ZIPI$kingdom<-as.factor(HR.BF.ZIPI$kingdom)
levels(HR.BF.ZIPI$kingdom)
HR.BF.ZIPI<-droplevels(HR.BF.ZIPI$kingdom)
#levels(HR.BF.ZIPI)<-c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#43A743","#97D1A0","grey","orange", "red","cyan","dodgerblue","black")
levels(HR.BF.ZIPI)

HR.BF.ZIPI<-read.csv("HR.Bac-Fun.ZIPI.csv",head=T,row.names = 1)
p2<-ggplot(HR.BF.ZIPI, aes(x = Pi, y = Zi, color=kingdom,shape=type)) +
  geom_point(aes(colour=kingdom,shape=type),size=5,alpha=0.5)+ 
  scale_color_manual(values = c("#ABDAEC","#43A743"))+
  scale_shape_manual(values = c(18,17,1))+
  labs(title="Cross-kingdom Network",subtitle="Symptomless Rhizosphere",x="Among-module connectivity (Pi)",y = "Within-module connectivity (Zi)")+
  scale_y_continuous(limits = c(-3,5))+
  scale_x_continuous(limits = c(0,1))+
  guides(fill=guide_legend(title= "kingdom"),shape=guide_legend(title= "type"))+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 20,face="bold"),
        legend.title = element_text(colour="black", size=20, face="bold"),
        legend.text = element_text(colour="black", size=16),
        axis.text.y=element_text(colour="black",size=16),
        axis.text.x=element_text(colour="black",size=16,angle = 0),
        axis.title=element_text(colour="black",size=20,face="bold"),
        plot.title=element_text(face = "bold", colour = "black", size = 20),
        plot.subtitle=element_text(colour = "black", size = 16))+
  geom_vline(xintercept = 0.62, color = 'gray34', linetype = 1, linewidth = 1.5) +
  geom_hline(yintercept = 2.5, color = 'gray34', linetype = 1, linewidth = 1.5) 
p2
ggsave("C:/graduate/WH赤霉病/pdf picture/cross-kingdom network HR ZIPI.pdf", p2,
       height=8,width =9,limitsize = FALSE )

rm(list=ls())
bac.coco<-read.csv("FZ.bac-fun.relabu.csv",head=T,row.names = 1)
matrix2igraph2<-function(matr,r.threshold,p.threshold){
  
  occor<-corAndPvalue(matr,method = c( "spearman")) 
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="fdr") 
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(matr)[2])
  occor.r<-occor$cor
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold]<-0 
  diag(occor.r)<-0  
  igraph<-graph.adjacency(occor.r,mode="undirected", weighted=TRUE, diag=FALSE)
  igraph<- igraph::simplify(igraph)
  igraph
}

bac<-matrix2igraph2(bac.coco,0.60,0.05)


bad.vs<-V(bac)[degree(bac) == 0] 
igraph <- delete.vertices(bac, bad.vs) 
igraph

##modularity 
E(igraph)$weight<-abs(E(igraph)$weight)
bac.fc<-cluster_fast_greedy(igraph)
modularity(bac.fc)##>0.4 indicate modular structures
modularity(bac.fc,membership(bac.fc))
membership(bac.fc)
sizes(bac.fc)


###Calculate vertex within-module degree z-score:???????????ġ?��?ӵ?
bac.tax_0.01<-read.csv("FZ.Bac-Fun.IDotutab.rar.csv",head=T,row.names=1)
bac.g.tax_0.01<-bac.tax_0.01[V(igraph)$name,]
bac.g.tax_0.01$phylum<-as.factor(bac.g.tax_0.01$phylum)
levels(bac.g.tax_0.01$phylum)
bac.g.tax<-droplevels(bac.g.tax_0.01$phylum)
bac.g.tax
levels(bac.g.tax)<-c("black","blue","dodgerblue","dodgerblue4","gray26",'red','orange')
  levels(bac.g.tax)

str(bac.g.tax_0.01)

library(brainGraph)

bac.comps <- membership(bac.fc)
bac_Pi<-part_coeff(igraph, membership(bac.fc))
bac.Zi<-within_module_deg_z_score(igraph, membership(bac.fc))
write.csv(bac_Pi,"FZ.bac-fun_Pi.csv")
write.csv(bac.Zi,"FZ.bac-fun_Zi.csv")

FZ.BF.ZIPI<-read.csv("FZ.Bac-Fun.ZIPI.csv",head=T,row.names = 1)
p3<-ggplot(FZ.BF.ZIPI, aes(x = Pi, y = Zi, color=kingdom,shape=type)) +
  geom_point(aes(colour=kingdom,shape=type),size=5,alpha=0.5)+ 
  scale_color_manual(values = c("#ABDAEC","#43A743"))+
  scale_shape_manual(values = c(18,17,1))+
  labs(title="Cross-kingdom Network",subtitle="Symptomatic Root",x="Among-module connectivity (Pi)",y = "Within-module connectivity (Zi)")+
  scale_y_continuous(limits = c(-3,5))+
  scale_x_continuous(limits = c(0,1))+
  guides(fill=guide_legend(title= "kingdom"),shape=guide_legend(title= "type"))+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 20,face="bold"),
        legend.title = element_text(colour="black", size=20, face="bold"),
        legend.text = element_text(colour="black", size=16),
        axis.text.y=element_text(colour="black",size=16),
        axis.text.x=element_text(colour="black",size=16,angle = 0),
        axis.title=element_text(colour="black",size=20,face="bold"),
        plot.title=element_text(face = "bold", colour = "black", size = 20),
        plot.subtitle=element_text(colour = "black", size = 16))+
  geom_vline(xintercept = 0.62, color = 'gray34', linetype = 1, linewidth = 1.5) +
  geom_hline(yintercept = 2.5, color = 'gray34', linetype = 1, linewidth = 1.5) 
p3
ggsave("C:/graduate/WH赤霉病/pdf picture/cross-kingdom network FZ ZIPI.pdf", p3,
       height=8,width =9,limitsize = FALSE )

rm(list=ls())
bac.coco<-read.csv("HZ.bac-fun.relabu.csv",head=T,row.names = 1)
matrix2igraph2<-function(matr,r.threshold,p.threshold){
  
  occor<-corAndPvalue(matr,method = c( "spearman")) 
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="fdr") 
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(matr)[2])
  occor.r<-occor$cor
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold]<-0 
  diag(occor.r)<-0  
  igraph<-graph.adjacency(occor.r,mode="undirected", weighted=TRUE, diag=FALSE)
  igraph<- igraph::simplify(igraph)
  igraph
}

bac<-matrix2igraph2(bac.coco,0.60,0.05)


bad.vs<-V(bac)[degree(bac) == 0] 
igraph <- delete.vertices(bac, bad.vs) 
igraph

##modularity 
E(igraph)$weight<-abs(E(igraph)$weight)
bac.fc<-cluster_fast_greedy(igraph)
modularity(bac.fc)##>0.4 indicate modular structures
modularity(bac.fc,membership(bac.fc))
membership(bac.fc)
sizes(bac.fc)


###Calculate vertex within-module degree z-score:???????????ġ?��?ӵ?
bac.tax_0.01<-read.csv("HZ.Bac-Fun.IDotutab.rar.csv",head=T,row.names=1)
bac.g.tax_0.01<-bac.tax_0.01[V(igraph)$name,]
bac.g.tax_0.01$phylum<-as.factor(bac.g.tax_0.01$phylum)
levels(bac.g.tax_0.01$phylum)
bac.g.tax<-droplevels(bac.g.tax_0.01$phylum)
bac.g.tax
levels(bac.g.tax)<-c("black","blue","dodgerblue","dodgerblue4","gray26",'red','orange')
  levels(bac.g.tax)

str(bac.g.tax_0.01)

library(brainGraph)

bac.comps <- membership(bac.fc)
bac_Pi<-part_coeff(igraph, membership(bac.fc))
bac.Zi<-within_module_deg_z_score(igraph, membership(bac.fc))
write.csv(bac_Pi,"HZ.bac-fun_Pi.csv")
write.csv(bac.Zi,"HZ.bac-fun_Zi.csv")

HZ.BF.ZIPI<-read.csv("HZ.Bac-Fun.ZIPI.csv",head=T,row.names = 1)
p4<-ggplot(HZ.BF.ZIPI, aes(x = Pi, y = Zi, color=kingdom,shape=type)) +
  geom_point(aes(colour=kingdom,shape=type),size=5,alpha=0.5)+ 
  scale_color_manual(values = c("#ABDAEC","#43A743"))+
  scale_shape_manual(values = c(18,17,1))+
  labs(title="Cross-kingdom Network",subtitle="Symptomless Root",x="Among-module connectivity (Pi)",y = "Within-module connectivity (Zi)")+
  scale_y_continuous(limits = c(-3,5))+
  scale_x_continuous(limits = c(0,1))+
  guides(fill=guide_legend(title= "kingdom"),shape=guide_legend(title= "type"))+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 20,face="bold"),
        legend.title = element_text(colour="black", size=20, face="bold"),
        legend.text = element_text(colour="black", size=16),
        axis.text.y=element_text(colour="black",size=16),
        axis.text.x=element_text(colour="black",size=16,angle = 0),
        axis.title=element_text(colour="black",size=20,face="bold"),
        plot.title=element_text(face = "bold", colour = "black", size = 20),
        plot.subtitle=element_text(colour = "black", size = 16))+
  geom_vline(xintercept = 0.62, color = 'gray34', linetype = 1, linewidth = 1.5) +
  geom_hline(yintercept = 2.5, color = 'gray34', linetype = 1, linewidth = 1.5) 
p4
ggsave("C:/graduate/WH赤霉病/pdf picture/cross-kingdom network HZ ZIPI.pdf", p4,
       height=8,width =9,limitsize = FALSE )

########################UpSet##########################
library(UpSetR)

setwd("C:/graduate/WH赤霉病/扩增子")
rm(list=ls())

otu1 <- read.csv("Rhizosphere Venn.csv",head=T,row.names = 1)
pdf("C:/graduate/WH赤霉病/pdf picture/Rhizospehre UpSet.pdf",height=8,width =12)
upset(otu1, nset = 6, nintersects = 100)
dev.off()

rm(list=ls())

otu2 <- read.csv("Root Venn.csv",head=T,row.names = 1)
pdf("C:/graduate/WH赤霉病/pdf picture/Root UpSet.pdf",height=8,width =12)
upset(otu2, nset = 6, nintersects = 100)
dev.off()
