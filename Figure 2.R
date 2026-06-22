################Bacterial network##############################
library(vegan)
library(psych)##
library(igraph)
library(beepr)

setwd("C:/graduate/WH赤霉病/扩增子/WH Bacteria")

rm(list=ls())
da.rar<-read.csv("da.rar.csv",head = T, row.names = 1)
env<-read.csv("env.csv",head = T, row.names = 1)
BB0<-read.csv("IDotutab.rar.csv", head = T, row.names = 1)
ID.tmp<-BB0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus", "OTU.ID")]

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
setwd("C:/graduate/WH赤霉病/扩增子/WH Bacteria")
BB0<-read.csv("IDotutab.rar.csv", head = T, row.names = 1)
ID.tmp<-BB0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus","OTU.ID")]
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
da.tmp$color[da.tmp$Cor>r.cutoff]<-"skyblue"
da.tmp$color[da.tmp$Cor< -r.cutoff]<-"grey"
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp1<-data.frame(da.tmp, M1, M2)

spman.r0 <- spman.d2
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[lower.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[lower.tri(Cor)]], Cor=Cor[lower.tri(Cor)])

P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[lower.tri(P0)]], 
                 col=colnames(P0)[col(P0)[lower.tri(P0)]], p=P0[lower.tri(P0)])

df <- data.frame(Cor.df,  P.df, Compartment = "Rhizosphere", State = "Symptomless")
da.tmp<-df.sig<- df[(df$Cor > r.cutoff | df$Cor < -r.cutoff) & df$p < p.cutoff,]
da.tmp$color[da.tmp$Cor>r.cutoff]<-"skyblue"
da.tmp$color[da.tmp$Cor< -r.cutoff]<-"grey"
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp2<-data.frame(da.tmp, M1, M2)


spman.r0 <- spman.d3
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[lower.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[lower.tri(Cor)]], Cor=Cor[lower.tri(Cor)])
P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[lower.tri(P0)]], 
                 col=colnames(P0)[col(P0)[lower.tri(P0)]], p=P0[lower.tri(P0)])

df <- data.frame(Cor.df,  P.df, Compartment = "Root", State = "Symptomatic")
da.tmp<-df.sig<- df[(df$Cor > r.cutoff | df$Cor < -r.cutoff) & df$p < p.cutoff,]
da.tmp$color[da.tmp$Cor>r.cutoff]<-"skyblue"
da.tmp$color[da.tmp$Cor< -r.cutoff]<-"grey" 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp3<-data.frame(da.tmp, M1, M2)

spman.r0 <- spman.d4
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[lower.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[lower.tri(Cor)]], Cor=Cor[lower.tri(Cor)])

P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[lower.tri(P0)]], 
                 col=colnames(P0)[col(P0)[lower.tri(P0)]], p=P0[lower.tri(P0)])
df <- data.frame(Cor.df,  P.df, Compartment = "Root", State = "Symptomless")
da.tmp<-df.sig<- df[(df$Cor > r.cutoff | df$Cor < -r.cutoff) & df$p < p.cutoff,]
da.tmp$color[da.tmp$Cor>r.cutoff]<-"skyblue"
da.tmp$color[da.tmp$Cor< -r.cutoff]<-"grey" 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp4<-data.frame(da.tmp, M1, M2)


####################
#Bacterial network##
####################

BB0<-read.csv("IDotutab.rar.csv", head = T, row.names = 1)
ID.tmp<-BB0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus","OTU.ID")]
ID.tmp$shape[ID.tmp$kingdom=="Bacteria"]<-"circle"; ID.tmp$shape[ID.tmp$kingdom=="Fungi"]<-"square"
ID.tmp$color <- "grey"
ID.tmp$color[ID.tmp$phylum=="Acidobacteriota"] <- "#807EBA"
ID.tmp$color[ID.tmp$phylum=="Actinobacteriota"] <- "#A7B7DF"
ID.tmp$color[ID.tmp$phylum=="Bacteroidota"] <- "#ABDAEC"
ID.tmp$color[ID.tmp$phylum=="Chloroflexi"] <- "#E9CEE5"
ID.tmp$color[ID.tmp$phylum=="Firmicutes"] <- "#FDD5C0"
ID.tmp$color[ID.tmp$phylum=="Gemmatimonadota"] <- "#FDD378"
ID.tmp$color[ID.tmp$phylum=="Myxococcota"] <- "lightcoral"
ID.tmp$color[ID.tmp$phylum=="Planctomycetota"] <- "rosybrown"
ID.tmp$color[ID.tmp$phylum=="Proteobacteria"] <- "#43A743"
ID.tmp$color[ID.tmp$phylum=="Verrucomicrobiota"] <- "#97D1A0"

da0<- df.tmp1
da0[is.na(da0)] = 0
da <- da0[(da0$kingdom == "Bacteria")& (da0$kingdom.1 == "Bacteria"),]

g <- graph.data.frame(da, directed=FALSE)
g.color = droplevels(ID.tmp[ID.tmp$OTU.ID %in% V(g)$name,])
g.color<-g.color[match(V(g)$name, g.color$OTU.ID),]
V(g)$color = as.character(g.color$color)
V(g)$shape <-as.character(g.color$shape)

num.edges = length(E(g)) 
num.vertices = length(V(g))
connectance = edge_density(g,loops=FALSE)# 
degree<-igraph::degree(g)
degree1<-as.data.frame(degree)
write.csv(degree1, "degeree BB Diseased Rhizosphere Network（正负）.csv")
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
df.tmp1<-data.frame(network = "BB", "Rhizosphere", "Symptomatic",  num.edges, num.vertices, connectance, average.degree, average.path.length, diameter, edge.connectivity, clustering.coefficient,
                    no.clusters, centralization.betweenness,centralization.degree,  Modularity, No.modules)
write.csv(df.tmp1, "BB Diseased Rhizosphere Network（正负）.csv")  


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

M1$BB<-da$color

E(g)$color = as.character(M1$BB)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/pdf picture/16S FR network.pdf",height=8,width =12)
plot(g, edge.width=0.5,  vertex.frame.color=NA,vertex.shape=c("circle"),vertex.label=NA,edge.lty=1, edge.curved=T,vertex.size=5,main="Symptomatic Rhizosphere",cex.main=0.5) 
legend(title="Bacterial Phylum",x=-1.95, y=0.8, 
       legend=c("Acidobacteriota","Actinobacteriota","Bacteroidota","Chloroflexi","Firmicutes","Gemmatimonadota","Myxococcota","Planctomycetota","Proteobacteria","Verrucomicrobiota","Other"), 
       pt.bg=c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#43A743","#97D1A0",  "grey"),
       col=c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#43A743","#97D1A0",  "grey"), 
       pch=21,pt.cex=1.5, cex=1.5, bty="n",box.lty=NA)
dev.off()

da0<- df.tmp2
da0[is.na(da0)] = 0
da <- da0[(da0$kingdom == "Bacteria")& (da0$kingdom.1 == "Bacteria"),]

g <- graph.data.frame(da, directed=FALSE)
g.color = droplevels(ID.tmp[ID.tmp$OTU.ID %in% V(g)$name,])
g.color<-g.color[match(V(g)$name, g.color$OTU.ID),]
V(g)$color = as.character(g.color$color)
V(g)$shape <-as.character(g.color$shape)

num.edges = length(E(g)) 
num.vertices = length(V(g))
connectance = edge_density(g,loops=FALSE)# 
degree<-igraph::degree(g)
degree2<-as.data.frame(degree)
write.csv(degree2, "degeree BB Healthy Rhizosphere Network（正负）.csv")
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
df.tmp2<-data.frame(network = "BB", "Rhizosphere", "Symptomless",  num.edges, num.vertices, connectance, average.degree, average.path.length, diameter, edge.connectivity, clustering.coefficient,
                    no.clusters, centralization.betweenness,centralization.degree,  Modularity, No.modules)
write.csv(df.tmp2, "BB Healthy Rhizosphere Network（正负）.csv")

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

M1$BB<-da$color

E(g)$color = as.character(M1$BB)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/pdf picture/16S HR network.pdf",height=8,width =12)
plot(g, edge.width=0.5,  vertex.frame.color=NA,vertex.shape=c("circle"),vertex.label=NA,edge.lty=1, edge.curved=T,vertex.size=5,main="Symptomless Rhizosphere",cex.main=0.5) 
legend(title="Bacterial Phylum",x=-1.95, y=0.8, 
       legend=c("Acidobacteriota","Actinobacteriota","Bacteroidota","Chloroflexi","Firmicutes","Gemmatimonadota","Myxococcota","Planctomycetota","Proteobacteria","Verrucomicrobiota","Other"), 
       pt.bg=c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#43A743","#97D1A0",  "grey"),
       col=c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#43A743","#97D1A0",  "grey"), 
       pch=21,pt.cex=1.5, cex=1.5, bty="n",box.lty=NA)
dev.off()

da0<- df.tmp3
da0[is.na(da0)] = 0
da <- da0[(da0$kingdom == "Bacteria")& (da0$kingdom.1 == "Bacteria"),]


g <- graph.data.frame(da, directed=FALSE)
g.color = droplevels(ID.tmp[ID.tmp$OTU.ID %in% V(g)$name,])
g.color<-g.color[match(V(g)$name, g.color$OTU.ID),]
V(g)$color = as.character(g.color$color)
V(g)$shape <-as.character(g.color$shape)



num.edges = length(E(g)) 
num.vertices = length(V(g))
connectance = edge_density(g,loops=FALSE)# 
degree<-igraph::degree(g)
degree3<-as.data.frame(degree)
write.csv(degree3, "degeree BB Diseased Root Network（正负）.csv")
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
df.tmp3<-data.frame(network = "BB", "Root", "Symptomatic",  num.edges, num.vertices, connectance, average.degree, average.path.length, diameter, edge.connectivity, clustering.coefficient,
                    no.clusters, centralization.betweenness,centralization.degree,  Modularity, No.modules)
write.csv(df.tmp3, "BB Diseased Root Network（正负）.csv")

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

M1$BB<-da$color

E(g)$color = as.character(M1$BB)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/pdf picture/16S FZ network.pdf",height=8,width =12)
plot(g, edge.width=0.5,  vertex.frame.color=NA,vertex.shape=c("circle"),vertex.label=NA,edge.lty=1, edge.curved=T,vertex.size=5,main="Symptomatic Root",cex.main=0.5) 
legend(title="Bacterial Phylum",x=-1.95, y=0.8, 
       legend=c("Acidobacteriota","Actinobacteriota","Bacteroidota","Chloroflexi","Firmicutes","Gemmatimonadota","Myxococcota","Planctomycetota","Proteobacteria","Verrucomicrobiota","Other"), 
       pt.bg=c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#43A743","#97D1A0",  "grey"),
       col=c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#43A743","#97D1A0",  "grey"), 
       pch=21,pt.cex=1.5, cex=1.5, bty="n",box.lty=NA)
dev.off()

da0<- df.tmp4
da0[is.na(da0)] = 0
da <- da0[(da0$kingdom == "Bacteria")& (da0$kingdom.1 == "Bacteria"),]

g <- graph.data.frame(da, directed=FALSE)
g.color = droplevels(ID.tmp[ID.tmp$OTU.ID %in% V(g)$name,])
g.color<-g.color[match(V(g)$name, g.color$OTU.ID),]
V(g)$color = as.character(g.color$color)
V(g)$shape <-as.character(g.color$shape)

num.edges = length(E(g)) 
num.vertices = length(V(g))
connectance = edge_density(g,loops=FALSE)# 
degree<-igraph::degree(g)
degree4<-as.data.frame(degree)
write.csv(degree4, "degeree BB Healthy Root Network（正负）.csv")
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
df.tmp4<-data.frame(network = "BB", "Root", "Symptomless",  num.edges, num.vertices, connectance, average.degree, average.path.length, diameter, edge.connectivity, clustering.coefficient,
                    no.clusters, centralization.betweenness,centralization.degree,  Modularity, No.modules)
write.csv(df.tmp4, "BB Healthy Root Network（正负）.csv")


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

M1$BB<-da$color
E(g)$color = as.character(M1$BB)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/pdf picture/16S HZ network.pdf",height=8,width =12)
plot(g, edge.width=0.5,  vertex.frame.color=NA,vertex.shape=c("circle"),vertex.label=NA,edge.lty=1, edge.curved=T,vertex.size=5,main="Symptomless Root",cex.main=0.5) 
legend(title="Bacterial Phylum",x=-1.95, y=0.8, 
       legend=c("Acidobacteriota","Actinobacteriota","Bacteroidota","Chloroflexi","Firmicutes","Gemmatimonadota","Myxococcota","Planctomycetota","Proteobacteria","Verrucomicrobiota","Other"), 
       pt.bg=c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#43A743","#97D1A0",  "grey"),
       col=c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378","lightcoral", "rosybrown", "#43A743","#97D1A0",  "grey"), 
       pch=21,pt.cex=1.5, cex=1.5, bty="n",box.lty=NA)
dev.off()

##############################Bacterial network properties########################
library(vegan)
library(ggplot2)
library(colorRamps)
library(ape)
library(splitstackshape)
library(reshape2)

setwd("C:/graduate/WH赤霉病/扩增子/WH Bacteria")
rm(list=ls())
BBNetwork <- read.csv("BB Network（正负）.csv", head = T, row.names = 1)

p<-ggplot(BBNetwork, aes(x = Compartment, y = num.edges, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=num.edges),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Bacterial Network", x = NULL, y = 'Number of edges', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 650))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S network num.edges.pdf", p,
       height=8,width =9,limitsize = FALSE )

p<-ggplot(BBNetwork, aes(x = Compartment, y = num.vertices, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=num.vertices),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Bacterial Network", x = NULL, y = 'Number of vertices', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 400))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S network num.vertices.pdf", p,
       height=8,width =9,limitsize = FALSE )

rm(list=ls())
degree <- read.csv("degree BB（正负）.csv",head=T)
degree.R<-degree[degree$Compartment=="Rhizosphere",] 
degree.Z<-degree[degree$Compartment=="Root",] 

shapiro.test(degree.R$degree)
shapiro.test(log(degree.R$degree))
shapiro.test(sqrt(degree.R$degree))
shapiro.test(degree.Z$degree)
shapiro.test(log(degree.Z$degree))
shapiro.test(sqrt(degree.Z$degree))

kruskal.test(degree.R$degree ~ degree.R$State)
kruskal.test(degree.Z$degree ~ degree.Z$State)

p<-ggplot(degree, aes(x = State, y = degree,color=State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.7,size=3.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "", x = NULL, y = 'Degree', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  scale_y_continuous(expand = c(0, 0), limit = c(0, 30))+
  ggtitle("Bacterial Network",expression('Rhizosphere:'~italic(χ) ^2~'= 2.757'^ns~~~~~'Root:'~italic(χ) ^2~'= 9.192**'))+
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
  facet_grid( .~Compartment)
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S network degree.pdf", p,
       height=8,width =9,limitsize = FALSE )


install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "impute","GO.db", "preprocessCore"))
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
install.packages(c("WGCNA", "stringr", "reshape2"), repos=site)

library(vegan)##decostand
library(WGCNA)##corAndPvalue
library(multtest)##mt.rawp2adjp
library(igraph)##graph.adjacency    delete.vertices

setwd("C:/graduate/WH赤霉病/扩增子/WH Bacteria")
rm(list=ls())


bac.coco<-read.csv("FR.bac.relabu.csv",head=T,row.names = 1)
#(bac.coco<-data.frame(decostand(bac.coco[,9:24],"hel")) 
#bac.coco<- t(bac.coco)


matrix2igraph2<-function(matr,r.threshold,p.threshold){
  
  occor<-corAndPvalue(matr,method = c( "spearman")) 
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="TSBH") 
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

bac.vs<-V(bac)[degree(bac) == 0] 
igraph <- delete.vertices(bac, bac.vs) 
igraph

a<-as.matrix(as_adjacency_matrix(igraph, type = c("both")))
write.csv(a,"FR.MATRIX.csv")

b<-read.csv("FR.MATRIX.csv",head=T,row.names=1)
network.raw<-b

sp.ra <-read.csv("FR.sp.ra.csv",head=F)
sp.ra<-as.matrix(sp.ra)
sp.ra2 <-read.csv("FR.sp.ra2.csv",head=F)
sp.ra2<-as.matrix(sp.ra2)


rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  
  remain.percent
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra2, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se,remains)
    names(result)<-c("remain.mean","remain.sd","remain.se","remains")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 Compartment=rep("Rhizosphere",40),State=rep("Symptomatic",4))

currentdat<-dat1

write.csv(currentdat,"16S.FR.rubustness100.csv")


rm(list=ls())
bac.coco<-read.csv("HR.bac.relabu.csv",head=T,row.names = 1)

matrix2igraph2<-function(matr,r.threshold,p.threshold){
  
  occor<-corAndPvalue(matr,method = c( "spearman")) 
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="TSBH") 
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


bac.vs<-V(bac)[degree(bac) == 0] 
igraph <- delete.vertices(bac, bac.vs) 
igraph

a<-as.matrix(as_adjacency_matrix(igraph, type = c("both")))
write.csv(a,"HR.MATRIX.csv")

b<-read.csv("HR.MATRIX.csv",head=T,row.names=1)
network.raw<-b

sp.ra <-read.csv("HR.sp.ra.csv",head=F)
sp.ra<-as.matrix(sp.ra)
sp.ra2 <-read.csv("HR.sp.ra2.csv",head=F)
sp.ra2<-as.matrix(sp.ra2)


rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  
  remain.percent
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra2, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se,remains)
    names(result)<-c("remain.mean","remain.sd","remain.se","remains")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 Compartment=rep("Rhizosphere",40),Statet=rep("Symptomless",4))

currentdat<-dat1

write.csv(currentdat,"16S.HR.rubustness100.csv")



rm(list=ls())
bac.coco<-read.csv("FZ.bac.relabu.csv",head=T,row.names = 1)

matrix2igraph2<-function(matr,r.threshold,p.threshold){
  
  occor<-corAndPvalue(matr,method = c( "spearman")) 
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="TSBH") 
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

bac.vs<-V(bac)[degree(bac) == 0]
igraph <- delete.vertices(bac, bac.vs) 
igraph

a<-as.matrix(as_adjacency_matrix(igraph, type = c("both")))
write.csv(a,"FZ.MATRIX.csv")

b<-read.csv("FZ.MATRIX.csv",head=T,row.names=1)
network.raw<-b

sp.ra <-read.csv("FZ.sp.ra.csv",head=F)
sp.ra<-as.matrix(sp.ra)
sp.ra2 <-read.csv("FZ.sp.ra2.csv",head=F)
sp.ra2<-as.matrix(sp.ra2)

rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  
  remain.percent
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra2, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se,remains)
    names(result)<-c("remain.mean","remain.sd","remain.se","remains")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 Compartment=rep("Root",40),State=rep("Symptomatic",4))

currentdat<-dat1

write.csv(currentdat,"16S.FZ.rubustness100.csv")



rm(list=ls())
bac.coco<-read.csv("HZ.bac.relabu.csv",head=T,row.names = 1)

matrix2igraph2<-function(matr,r.threshold,p.threshold){
  
  occor<-corAndPvalue(matr,method = c( "spearman")) 
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="TSBH") 
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

bac.vs<-V(bac)[degree(bac) == 0] 
igraph <- delete.vertices(bac, bac.vs) 
igraph

a<-as.matrix(as_adjacency_matrix(igraph, type = c("both")))
write.csv(a,"HZ.MATRIX.csv")

b<-read.csv("HZ.MATRIX.csv",head=T,row.names=1)
network.raw<-b


sp.ra <-read.csv("HZ.sp.ra.csv",head=F)
sp.ra<-as.matrix(sp.ra)
sp.ra2 <-read.csv("HZ.sp.ra2.csv",head=F)
sp.ra2<-as.matrix(sp.ra2)

rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  
  remain.percent
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra2, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se,remains)
    names(result)<-c("remain.mean","remain.sd","remain.se","remains")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 Compartment=rep("Root",40),State=rep("Symptomless",4))

currentdat<-dat1

write.csv(currentdat,"16S.HZ.rubustness100.csv")

library(ggplot2)
All <- read.csv("16S.All.ANNOVA.ROBUST.csv",head=T)
Rhizosphere<-All[All$Compartment=="Rhizosphere",]
Root<-All[All$Compartment=="Root",]


shapiro.test(Rhizosphere$Robust) 
shapiro.test(log(Rhizosphere$Robust)) 
shapiro.test(sqrt(Rhizosphere$Robust)) 

bartlett.test(Robust~ interaction(State),data=Rhizosphere)

anova_Rhizosphere <- aov(Robust ~ State,data = Rhizosphere)
summary(anova_Rhizosphere)


shapiro.test(Root$Robust) 
shapiro.test(log(Root$Robust)) 
shapiro.test(sqrt(Root$Robust)) 


bartlett.test(Robust~ interaction(State),data=Root) 

anova_Root <- aov(Robust ~ State,data = Root)
summary(anova_Root)

p<-ggplot(All, aes(x = State, y = Robust,color=State)) +
  geom_boxplot(size=0.8, width = 0.3,alpha=0) +
  geom_jitter(position =position_jitter(0.2),alpha=0.2,size=1.5)+
  scale_color_manual(values=c("#8b66b8","#51c4c2"))+
  labs(title = "", x = NULL, y = 'Robustness', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  scale_y_continuous(expand = c(0, 0), limit = c(0.4, 0.5))+
  ggtitle("Bacterial Network",expression('Rhizosphere:'~italic(F)~'= 53.700***'~~~~~'Root:'~italic(F)~'= 204.000***'))+
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
  facet_grid( .~All$Compartment)
p
ggsave("C:/graduate/WH赤霉病/pdf picture/16S network Robustness.pdf", p,
       height=8,width =9,limitsize = FALSE )

######################Fungal network#####################
library(vegan)
library(psych)##
library(igraph)
library(beepr)

setwd("C:/graduate/WH赤霉病/扩增子/WH Fungi")
rm(list=ls())

da.rar<-read.csv("da.rar.csv",head = T, row.names = 1)
env<-read.csv("env.csv",head = T, row.names = 1)
FF0<-read.csv("IDotutab.rar.csv", head = T, row.names = 1)
ID.tmp<-FF0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus", "OTU.ID")]

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
setwd("C:/graduate/WH赤霉病/扩增子/WH Fungi")
FF0<-read.csv("IDotutab.rar.csv", head = T, row.names = 1)
ID.tmp<-FF0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus","OTU.ID")]
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
da.tmp$color[da.tmp$Cor>r.cutoff]<-"red"
da.tmp$color[da.tmp$Cor< -r.cutoff]<-"grey" 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp1<-data.frame(da.tmp, M1, M2)


spman.r0 <- spman.d2
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[lower.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[lower.tri(Cor)]], Cor=Cor[lower.tri(Cor)])

P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[lower.tri(P0)]], 
                 col=colnames(P0)[col(P0)[lower.tri(P0)]], p=P0[lower.tri(P0)])

df <- data.frame(Cor.df,  P.df, Compartment = "Rhizosphere", State = "Symptomless")
da.tmp<-df.sig<- df[ (df$Cor > r.cutoff | df$Cor < -r.cutoff) & df$p < p.cutoff,] 
da.tmp$color[da.tmp$Cor>r.cutoff]<-"red"
da.tmp$color[da.tmp$Cor< -r.cutoff]<-"grey" 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp2<-data.frame(da.tmp, M1, M2)


spman.r0 <- spman.d3
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[lower.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[lower.tri(Cor)]], Cor=Cor[lower.tri(Cor)])
P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[lower.tri(P0)]], 
                 col=colnames(P0)[col(P0)[lower.tri(P0)]], p=P0[lower.tri(P0)])

df <- data.frame(Cor.df,  P.df, Compartment = "Root", State = "Symptomatic")
da.tmp<-df.sig<- df[ (df$Cor > r.cutoff | df$Cor < -r.cutoff) & df$p < p.cutoff,] 
da.tmp$color[da.tmp$Cor>r.cutoff]<-"red"
da.tmp$color[da.tmp$Cor< -r.cutoff]<-"grey" 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp3<-data.frame(da.tmp, M1, M2)

spman.r0 <- spman.d4
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[lower.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[lower.tri(Cor)]], Cor=Cor[lower.tri(Cor)])

P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[lower.tri(P0)]], 
                 col=colnames(P0)[col(P0)[lower.tri(P0)]], p=P0[lower.tri(P0)])
df <- data.frame(Cor.df,  P.df, Compartment = "Root", State = "Symptomless")
da.tmp<-df.sig<- df[ (df$Cor > r.cutoff | df$Cor < -r.cutoff) & df$p < p.cutoff,] 
da.tmp$color[da.tmp$Cor>r.cutoff]<-"red"
da.tmp$color[da.tmp$Cor< -r.cutoff]<-"grey" 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp4<-data.frame(da.tmp, M1, M2)


####################
#Fungal network##
####################
###Phylum###
FF0<-read.csv("IDotutab.rar.csv", head = T, row.names = 1)
ID.tmp<-FF0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus","OTU.ID")]
ID.tmp$shape[ID.tmp$kingdom=="Bacteria"]<-"circle"; ID.tmp$shape[ID.tmp$kingdom=="Fungi"]<-"square"
ID.tmp$color <- "grey"
ID.tmp$color[ID.tmp$phylum=="Ascomycota"] <- "#807EBA"
ID.tmp$color[ID.tmp$phylum=="Basidiomycota"] <- "#ABDAEC"
ID.tmp$color[ID.tmp$phylum=="Chytridiomycota"] <- "#E9CEE5"
ID.tmp$color[ID.tmp$phylum=="Entomophthoromycota"] <- "#A7B7DF"
ID.tmp$color[ID.tmp$phylum=="Kickxellomycota"] <- "#FDD5C0"
ID.tmp$color[ID.tmp$phylum=="Mortierellomycota"] <- "#FDD378"
ID.tmp$color[ID.tmp$phylum=="Mucoromycota"] <- "lightcoral"
ID.tmp$color[ID.tmp$phylum=="Olpidiomycota"] <- "rosybrown"
ID.tmp$color[ID.tmp$phylum=="Rozellomycota"] <- "#43A743"
ID.tmp$color[ID.tmp$phylum=="Zoopagomycota"] <- "#97D1A0"


da0<- df.tmp1
da0[is.na(da0)] = 0
da <- da0[da0$kingdom == "Fungi"&da0$kingdom.1 == "Fungi",]

g <- graph.data.frame(da, directed=FALSE)
g.color = droplevels(ID.tmp[ID.tmp$OTU.ID %in% V(g)$name,])
g.color<-g.color[match(V(g)$name, g.color$OTU.ID),]
V(g)$color = as.character(g.color$color)
V(g)$shape <-as.character(g.color$shape)

num.edges = length(E(g)) 
num.vertices = length(V(g))
connectance = edge_density(g,loops=FALSE)# 
degree<-igraph::degree(g)
degree1<-as.data.frame(degree)
write.csv(degree1, "degeree FF Diseased Rhizosphere Network（正负）.csv")
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
df.tmp1<-data.frame(network = "FF", "Rhizosphere", "Symptomatic",  num.edges, num.vertices, connectance, average.degree, average.path.length, diameter, edge.connectivity, clustering.coefficient,
                    no.clusters, centralization.betweenness,centralization.degree,  Modularity, No.modules)
write.csv(df.tmp1, "FF Diseased Rhizosphere Network（正负）.csv")

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

M1$FF<-da$color

E(g)$color = as.character(M1$FF)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/pdf picture/ITS FR network.pdf",height=8,width =12)
plot(g, edge.width=0.3,  vertex.frame.color=NA,vertex.shape=c("square"),vertex.label=NA,edge.lty=1, edge.curved=T,vertex.size=5,main="Symptomatic Rhizosphere",cex.main=0.5) 
legend(title="Fungal Phylum",x=-1.95, y=0.5, 
       title.adj = 0.4,
       legend=c("Ascomycota","Basidiomycota","Chytridiomycota","Mortierellomycota","Mucoromycota","Olpidiomycota","Other"), 
       pt.bg=c("#807EBA","#ABDAEC","#E9CEE5","#FDD378","lightcoral", "rosybrown",  "grey"),
       col=c("#807EBA","#ABDAEC","#E9CEE5","#FDD378","lightcoral", "rosybrown",  "grey"), 
       pch=22,pt.cex=1.5,cex=1.5, bty="n",box.lty=NA)
dev.off()

da0<- df.tmp2
da0[is.na(da0)] = 0
da <- da0[da0$kingdom == "Fungi"&da0$kingdom.1 == "Fungi",]

g <- graph.data.frame(da, directed=FALSE)
g.color = droplevels(ID.tmp[ID.tmp$OTU.ID %in% V(g)$name,])
g.color<-g.color[match(V(g)$name, g.color$OTU.ID),]
V(g)$color = as.character(g.color$color)
V(g)$shape <-as.character(g.color$shape)

num.edges = length(E(g)) 
num.vertices = length(V(g))
connectance = edge_density(g,loops=FALSE)# 
degree<-igraph::degree(g)
degree2<-as.data.frame(degree)
write.csv(degree2, "degeree FF Healthy Rhizosphere Network（正负）.csv")
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
df.tmp2<-data.frame(network = "FF", "Rhizosphere", "Symptomless",  num.edges, num.vertices, connectance, average.degree, average.path.length, diameter, edge.connectivity, clustering.coefficient,
                    no.clusters, centralization.betweenness,centralization.degree,  Modularity, No.modules)
write.csv(df.tmp2, "FF Healthy Rhizosphere Network（正负）.csv")

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

M1$FF<-da$color

E(g)$color = as.character(M1$FF)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/pdf picture/ITS HR network.pdf",height=8,width =12)
plot(g, edge.width=0.2,  vertex.frame.color=NA,vertex.shape=c("square"),vertex.label=NA,edge.lty=1, edge.curved=T,vertex.size=5,main="Symptomless Rhizosphere",cex.main=0.5) 
legend(title="Fungal Phylum",x=-1.95, y=0.5, 
       title.adj = 0.4,
       legend=c("Ascomycota","Basidiomycota","Chytridiomycota","Mortierellomycota","Mucoromycota","Olpidiomycota","Other"), 
       pt.bg=c("#807EBA","#ABDAEC","#E9CEE5","#FDD378","lightcoral", "rosybrown",  "grey"),
       col=c("#807EBA","#ABDAEC","#E9CEE5","#FDD378","lightcoral", "rosybrown",  "grey"), 
       pch=22,pt.cex=1.5,cex=1.5, bty="n",box.lty=NA)
dev.off()

da0<- df.tmp3
da0[is.na(da0)] = 0
da <- da0[da0$kingdom == "Fungi"&da0$kingdom.1 == "Fungi",]

g <- graph.data.frame(da, directed=FALSE)
g.color = droplevels(ID.tmp[ID.tmp$OTU.ID %in% V(g)$name,])
g.color<-g.color[match(V(g)$name, g.color$OTU.ID),]
V(g)$color = as.character(g.color$color)
V(g)$shape <-as.character(g.color$shape)


num.edges = length(E(g)) 
num.vertices = length(V(g))
connectance = edge_density(g,loops=FALSE)# 
degree<-igraph::degree(g)
degree3<-as.data.frame(degree)
write.csv(degree3, "degeree FF Diseased Root Network（正负）.csv")
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
df.tmp3<-data.frame(network = "FF", "Root", "Symptomatic",  num.edges, num.vertices, connectance, average.degree, average.path.length, diameter, edge.connectivity, clustering.coefficient,
                    no.clusters, centralization.betweenness,centralization.degree,  Modularity, No.modules)
write.csv(df.tmp3, "FF Diseased Root Network（正负）.csv")


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

M1$FF<-da$color

E(g)$color = as.character(M1$FF)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/pdf picture/ITS FZ network.pdf",height=8,width =12)
plot(g, edge.width=0.3,  vertex.frame.color=NA,vertex.shape=c("square"),vertex.label=NA,edge.lty=1, edge.curved=T,vertex.size=5,main="Symptomatic Root",cex.main=0.5) 
legend(title="Fungal Phylum",x=-1.95, y=0.5, 
       title.adj = 0.4,
       legend=c("Ascomycota","Basidiomycota","Chytridiomycota","Mortierellomycota","Mucoromycota","Olpidiomycota","Other"), 
       pt.bg=c("#807EBA","#ABDAEC","#E9CEE5","#FDD378","lightcoral", "rosybrown",  "grey"),
       col=c("#807EBA","#ABDAEC","#E9CEE5","#FDD378","lightcoral", "rosybrown",  "grey"), 
       pch=22,pt.cex=1.5,cex=1.5, bty="n",box.lty=NA)
dev.off()

da0<- df.tmp4
da0[is.na(da0)] = 0
da <- da0[da0$kingdom == "Fungi"&da0$kingdom.1 == "Fungi",]

g <- graph.data.frame(da, directed=FALSE)
g.color = droplevels(ID.tmp[ID.tmp$OTU.ID %in% V(g)$name,])
g.color<-g.color[match(V(g)$name, g.color$OTU.ID),]
V(g)$color = as.character(g.color$color)
V(g)$shape <-as.character(g.color$shape)

num.edges = length(E(g)) 
num.vertices = length(V(g))
connectance = edge_density(g,loops=FALSE)# 
degree<-igraph::degree(g)
degree4<-as.data.frame(degree)
write.csv(degree4, "degeree FF Healthy Root Network（正负）.csv")
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
df.tmp4<-data.frame(network = "FF", "Root", "Symptomless",  num.edges, num.vertices, connectance, average.degree, average.path.length, diameter, edge.connectivity, clustering.coefficient,
                    no.clusters, centralization.betweenness,centralization.degree,  Modularity, No.modules)
write.csv(df.tmp4, "FF Healthy Root Network（正负）.csv")


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

M1$FF<-da$color
E(g)$color = as.character(M1$FF)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/pdf picture/ITS HZ network.pdf",height=8,width =12)
plot(g, edge.width=0.3,  vertex.frame.color=NA,vertex.shape=c("square"),vertex.label=NA,edge.lty=1, edge.curved=T,vertex.size=5,main="Symptomless Root",cex.main=0.5) 
legend(title="Fungal Phylum",x=-1.95, y=0.5, 
       title.adj = 0.4,
       legend=c("Ascomycota","Basidiomycota","Chytridiomycota","Mortierellomycota","Mucoromycota","Olpidiomycota","Other"), 
       pt.bg=c("#807EBA","#ABDAEC","#E9CEE5","#FDD378","lightcoral", "rosybrown",  "grey"),
       col=c("#807EBA","#ABDAEC","#E9CEE5","#FDD378","lightcoral", "rosybrown",  "grey"), 
       pch=22,pt.cex=1.5,cex=1.5, bty="n",box.lty=NA)
dev.off()

##########################Fungal network properties#########################
library(vegan)
library(ggplot2)
library(colorRamps)
library(ape)
library(splitstackshape)
library(reshape2)

setwd("C:/graduate/WH赤霉病/扩增子/WH Fungi")

rm(list=ls())
FFNetwork <- read.csv("FF Network（正负）.csv", head = T, row.names = 1)

p<-ggplot(FFNetwork, aes(x = Compartment, y = num.edges, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=num.edges),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Fungal Network", x = NULL, y = 'Number of edges', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 250))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS network num.edges.pdf", p,
       height=8,width =9,limitsize = FALSE )

p<-ggplot(FFNetwork, aes(x = Compartment, y = num.vertices, fill = State)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.5,alpha = 0.5) +
  #geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, size = 0.3, position = position_dodge(0.9))# +
  geom_text(aes(label=num.vertices),position=position_dodge(width = 0.9),size=5,vjust=-0.25 ) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "Fungal Network", x = NULL, y = 'Number of vertices', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.position = c(0.9, 0.85)) +
  theme_bw()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 16),
        axis.text.y =element_text(color="black",size = 16),
        #legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(colour="black", size=16),
        legend.title = element_text(face = "bold", colour = "black", size = 20),
        plot.title=element_text(face = "bold", colour = "black", size = 20)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 150))
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS network num.vertices.pdf", p,
       height=8,width =9,limitsize = FALSE )

rm(list=ls())
degree <- read.csv("degree FF（正负）.csv",head=T)
degree.R<-degree[degree$Compartment=="Rhizosphere",] 
degree.Z<-degree[degree$Compartment=="Root",] 

shapiro.test(degree.R$degree) #no
shapiro.test(log(degree.R$degree)) #no
shapiro.test(sqrt(degree.R$degree))# no
shapiro.test(degree.Z$degree) #no
shapiro.test(log(degree.Z$degree)) #no
shapiro.test(sqrt(degree.Z$degree))#no

kruskal.test(degree.R$degree ~ degree.R$State)
kruskal.test(degree.Z$degree ~ degree.Z$State)

p<-ggplot(degree, aes(x = State, y = degree,color=State)) +
  geom_violin(position = position_dodge(width = 0.3),alpha=0.2,color="white") +
  geom_jitter(aes(color = State),position =position_jitter(0.3),alpha=0.7,size=3.5) +
  scale_color_manual(values = c("#8b66b8","#51c4c2")) +
  scale_fill_manual(values = c("#8b66b8","#51c4c2")) +
  labs(title = "", x = NULL, y = 'Degree', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  scale_y_continuous(expand = c(0, 0), limit = c(0, 20))+
  ggtitle("Fungal Network",expression('Rhizosphere:'~italic(χ) ^2~'= 8.504**'~~~~~'Root:'~italic(χ) ^2~'= 7.980**'))+
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
  facet_grid( .~Compartment)
p
ggsave("C:/graduate/WH赤霉病/pdf picture/ITS network degree.pdf", p,
       height=8,width =9,limitsize = FALSE )


install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "impute","GO.db", "preprocessCore"))
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
install.packages(c("WGCNA", "stringr", "reshape2"), repos=site)

library(vegan)##decostand
library(WGCNA)##corAndPvalue
library(multtest)##mt.rawp2adjp
library(igraph)##graph.adjacency    delete.vertices确

setwd("C:/graduate/WH赤霉病/扩增子/WH Fungi")
rm(list = ls())

sp.ra <-read.csv("FR.sp.ra.csv",head=F)
sp.ra<-as.matrix(sp.ra)
sp.ra2 <-read.csv("FR.sp.ra2.csv",head=F)
sp.ra2<-as.matrix(sp.ra2)

fun.coco<-read.csv("FR.fung.relabu.csv",head=T,row.names = 1)

matrix2igraph2<-function(matr,r.threshold,p.threshold){
  
  occor<-corAndPvalue(matr,method = c( "spearman")) 
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="TSBH") 
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(matr)[2])
  occor.r<-occor$cor
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold]<-0 
  diag(occor.r)<-0  
  igraph<-graph.adjacency(occor.r,mode="undirected", weighted=TRUE, diag=FALSE)
  igraph<- igraph::simplify(igraph)
  igraph
}

fun<-matrix2igraph2(fun.coco,0.60,0.05) 


fun.vs<-V(fun)[degree(fun) == 0] 
igraph <- delete.vertices(fun, fun.vs) 
igraph

a<-as.matrix(as_adjacency_matrix(igraph, type = c("both")))
write.csv(a,"FR.MATRIX.csv")

b<-read.csv("FR.MATRIX.csv",head=T,row.names=1)
network.raw<-b


rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  
  remain.percent
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra2, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se,remains)
    names(result)<-c("remain.mean","remain.sd","remain.se","remains")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 Compartment=rep("Rhizosphere",40),State=rep("Symptomatic",4))

currentdat<-dat1

write.csv(currentdat,"ITS.FR.rubustness100.csv")


sp.ra <-read.csv("HR.sp.ra.csv",head=F)
sp.ra<-as.matrix(sp.ra)
sp.ra2 <-read.csv("HR.sp.ra2.csv",head=F)
sp.ra2<-as.matrix(sp.ra2)

fun.coco<-read.csv("HR.fung.relabu.csv",head=T,row.names = 1)

matrix2igraph2<-function(matr,r.threshold,p.threshold){
  
  occor<-corAndPvalue(matr,method = c( "spearman")) 
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="TSBH") 
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(matr)[2])
  occor.r<-occor$cor
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold]<-0 
  diag(occor.r)<-0  
  igraph<-graph.adjacency(occor.r,mode="undirected", weighted=TRUE, diag=FALSE)
  igraph<- igraph::simplify(igraph)
  igraph
}

fun<-matrix2igraph2(fun.coco,0.60,0.05) 

fun.vs<-V(fun)[degree(fun) == 0] 
igraph <- delete.vertices(fun, fun.vs) 
igraph

a<-as.matrix(as_adjacency_matrix(igraph, type = c("both")))
write.csv(a,"HR.MATRIX.csv")

b<-read.csv("HR.MATRIX.csv",head=T,row.names=1)
network.raw<-b


rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  
  remain.percent
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra2, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se,remains)
    names(result)<-c("remain.mean","remain.sd","remain.se","remains")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 Compartment=rep("Rhizosphere",40),State=rep("Symptomless",4))

currentdat<-dat1

write.csv(currentdat,"ITS.HR.rubustness100.csv")


sp.ra <-read.csv("FZ.sp.ra.csv",head=F)
sp.ra<-as.matrix(sp.ra)
sp.ra2 <-read.csv("FZ.sp.ra2.csv",head=F)
sp.ra2<-as.matrix(sp.ra2)

fun.coco<-read.csv("FZ.fung.relabu.csv",head=T,row.names = 1)

matrix2igraph2<-function(matr,r.threshold,p.threshold){
  
  occor<-corAndPvalue(matr,method = c( "spearman")) 
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="TSBH") 
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(matr)[2])
  occor.r<-occor$cor
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold]<-0 
  diag(occor.r)<-0  
  igraph<-graph.adjacency(occor.r,mode="undirected", weighted=TRUE, diag=FALSE)
  igraph
}

fun<-matrix2igraph2(fun.coco,0.60,0.05)

fun.vs<-V(fun)[degree(fun) == 0] 
igraph <- delete.vertices(fun, fun.vs) 
igraph

a<-as.matrix(as_adjacency_matrix(igraph, type = c("both")))
write.csv(a,"FZ.MATRIX.csv")

b<-read.csv("FZ.MATRIX.csv",head=T,row.names=1)
network.raw<-b

rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  
  remain.percent
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra2, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se,remains)
    names(result)<-c("remain.mean","remain.sd","remain.se","remains")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 Compartment=rep("Root",40),State=rep("Symptomatic",4))

currentdat<-dat1

write.csv(currentdat,"ITS.FZ.rubustness100.csv")


sp.ra <-read.csv("HZ.sp.ra.csv",head=F)
sp.ra<-as.matrix(sp.ra)
sp.ra2 <-read.csv("HZ.sp.ra2.csv",head=F)
sp.ra2<-as.matrix(sp.ra2)

fun.coco<-read.csv("HZ.fung.relabu.csv",head=T,row.names = 1)

matrix2igraph2<-function(matr,r.threshold,p.threshold){
  
  occor<-corAndPvalue(matr,method = c( "spearman")) 
  mtadj<-mt.rawp2adjp(unlist(occor$p),proc="TSBH") 
  adpcor<-mtadj$adjp[order(mtadj$index),2]
  occor.p<-matrix(adpcor,dim(matr)[2])
  occor.r<-occor$cor
  occor.r[occor.p>p.threshold|abs(occor.r)<r.threshold]<-0 
  diag(occor.r)<-0  
  igraph<-graph.adjacency(occor.r,mode="undirected", weighted=TRUE, diag=FALSE)
  igraph<- igraph::simplify(igraph)
  igraph
}

fun<-matrix2igraph2(fun.coco,0.60,0.05)


fun.vs<-V(fun)[degree(fun) == 0] 
igraph <- delete.vertices(fun, fun.vs) 
igraph

a<-as.matrix(as_adjacency_matrix(igraph, type = c("both")))
write.csv(a,"HZ.MATRIX.csv")

b<-read.csv("HZ.MATRIX.csv",head=T,row.names=1)
network.raw<-b


rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  
  remain.percent
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra2, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se,remains)
    names(result)<-c("remain.mean","remain.sd","remain.se","remains")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 Compartment=rep("Root",40),State=rep("Symptomless",4))

currentdat<-dat1

write.csv(currentdat,"ITS.HZ.rubustness100.csv")

library(ggplot2)
All <- read.csv("ITS.All.ANNOVA.ROBUST.csv",head=T)
Rhizosphere<-All[All$Compartment=="Rhizosphere",]
Root<-All[All$Compartment=="Root",]

shapiro.test(Rhizosphere$Robust) 
shapiro.test(log(Rhizosphere$Robust)) 
shapiro.test(sqrt(Rhizosphere$Robust)) 

kruskal.test(Rhizosphere$Robust ~ Rhizosphere$State)


shapiro.test(Root$Robust) 
shapiro.test(log(Root$Robust)) 
shapiro.test(sqrt(Root$Robust)) 

bartlett.test(Robust~ interaction(State),data=Root)

anova_Root <- aov(Robust ~ State,data = Root)
summary(anova_Root)



p<-ggplot(All, aes(x = State, y = Robust,color=State)) +
  geom_boxplot(size=0.8, width = 0.3,alpha=0) +
  geom_jitter(position =position_jitter(0.2),alpha=0.2,size=1.5)+
  scale_color_manual(values=c("#8b66b8","#51c4c2"))+
  labs(title = "", x = NULL, y = 'Robustness', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'))+
  scale_y_continuous(expand = c(0, 0), limit = c(0.44, 0.5))+
  ggtitle("Fungal Network",expression('Rhizosphere:'~italic(蠂) ^2~'= 141.500***'~~~~~'Root:'~italic(F)~'= 5.000*'))+
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
  facet_grid( .~All$Compartment)
p
ggsave("C:/graduate/WH璧ら湁鐥?/pdf picture/ITS network Robustness.pdf", p,
       height=8,width =9,limitsize = FALSE )
