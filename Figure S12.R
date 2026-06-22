################Bacterial network module#################
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
# Module bac-bac###
####################   
BB0<-read.csv("IDotutab.rar.csv", head = T, row.names = 1)
ID.tmp<-BB0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus",	"OTU.ID")]

da0<- df.tmp1
da0[is.na(da0)] = 0
da <- da0[(da0$kingdom == "Bacteria")& (da0$kingdom.1 == "Bacteria"),]
g <- graph.data.frame(da, directed=FALSE)

M1$BB<-da$color
E(g)$color = as.character(M1$BB)

fun.fc<-cluster_fast_greedy(g)
print(modularity(fun.fc))
modularity(fun.fc,membership(fun.fc))
membership(fun.fc)
print(sizes(fun.fc))###
fun.comps <- membership(fun.fc)
colbar <-c("#0000FF","#FF3399","#FFCC33","#9999FF", "#FF99FF", "#FFFFCC",
           "#33FFFF","#CCFF00","#CCFFCC","#FF3300","#CC99CC","#CC66CC","#FFFFCC","#CC0033","#666666","slateblue1","springgreen","steelblue1","tan1","thistle1","tomato","turquoise", "violet", "red", "yellowgreen",    "peachpuff", "peru", "pink", "plum2", "purple","wheat", "cornsilk3","cornsilk","coral4","coral",
           "chocolate4","chocolate","black","chartreuse4","chartreuse","burlywood4","burlywood","brown4","blue","bisque4","bisque",
           "azure4","azure","aquamarine4","aquamarine","antiquewhite4","antiquewhite","aliceblue","dodgerblue4","dodgerblue","dimgrey",
           "deepskyblue4", "deepskyblue", "deeppink4", "deeppink","darkviolet", "darkslategray4","darkslategray",
           "darkseagreen4", "darkseagreen", "darksalmon", "darkred", "darkorchid4", "darkorchid", "darkorange4",
           "darkorange","firebrick","darkgreen", "darkgoldenrod4","grey", "darkgoldenrod", "darkcyan","orange","rosybrown2","thistle1","lavender","paleturquoise","darkseagreen2","mistyrose","slateblue","darkolivegreen1","palevioletred1")
V(g)$color <- colbar[fun.comps]

set.seed(123)
pdf("C:/graduate/WH赤霉病/pdf picture/16S FR network module.pdf",height=8,width =9)
plot(g, edge.width=0.5, vertex.frame.color=NA,vertex.label=NA,edge.lty=1, edge.curved=T,vertex.size=3,margin=c(0, 0,0,0),main="Symptomatic Rhizosphere",cex.main=0.5) 
dev.off()

da0<- df.tmp2
da0[is.na(da0)] = 0
da <- da0[(da0$kingdom == "Bacteria")& (da0$kingdom.1 == "Bacteria"),]
g <- graph.data.frame(da, directed=FALSE)

M1$BB<-da$color
E(g)$color = as.character(M1$BB)

fun.fc<-cluster_fast_greedy(g)
print(modularity(fun.fc))
modularity(fun.fc,membership(fun.fc))
membership(fun.fc)
print(sizes(fun.fc))###
fun.comps <- membership(fun.fc)
colbar <-c("#0000FF","#FF3399","#FFCC33","#9999FF", "#FF99FF", "#FFFFCC",
           "#33FFFF","#CCFF00","#CCFFCC","#FF3300","#CC99CC","#CC66CC","#FFFFCC","#CC0033","#666666","slateblue1","springgreen","steelblue1","tan1","thistle1","tomato","turquoise", "violet", "red", "yellowgreen",    "peachpuff", "peru", "pink", "plum2", "purple","wheat", "cornsilk3","cornsilk","coral4","coral",
           "chocolate4","chocolate","black","chartreuse4","chartreuse","burlywood4","burlywood","brown4","blue","bisque4","bisque",
           "azure4","azure","aquamarine4","aquamarine","antiquewhite4","antiquewhite","aliceblue","dodgerblue4","dodgerblue","dimgrey",
           "deepskyblue4", "deepskyblue", "deeppink4", "deeppink","darkviolet", "darkslategray4","darkslategray",
           "darkseagreen4", "darkseagreen", "darksalmon", "darkred", "darkorchid4", "darkorchid", "darkorange4",
           "darkorange","firebrick","darkgreen", "darkgoldenrod4","grey", "darkgoldenrod", "darkcyan")
V(g)$color <- colbar[fun.comps]

set.seed(123)
pdf("C:/graduate/WH赤霉病/pdf picture/16S HR network module.pdf",height=8,width =9)
plot(g, edge.width=0.5, vertex.frame.color=NA,vertex.label=NA,edge.lty=1,edge.curved=T,vertex.size=3,margin=c(0, 0,0,0),main="Symptomless Rhizosphere",cex.main=0.5) 
dev.off()

da0<- df.tmp3
da0[is.na(da0)] = 0
da <- da0[(da0$kingdom == "Bacteria")& (da0$kingdom.1 == "Bacteria"),]
g <- graph.data.frame(da, directed=FALSE)

M1$BB<-da$color
E(g)$color = as.character(M1$BB)

fun.fc<-cluster_fast_greedy(g)
print(modularity(fun.fc))
modularity(fun.fc,membership(fun.fc))
membership(fun.fc)
print(sizes(fun.fc))###
fun.comps <- membership(fun.fc)
colbar <-c("#0000FF","#FF3399","#FFCC33","#9999FF", "#FF99FF", "#FFFFCC",
           "#33FFFF","#CCFF00","#CCFFCC","#FF3300","#CC99CC","#CC66CC","#FFFFCC","#CC0033","#666666","slateblue1","springgreen","steelblue1","tan1","thistle1","tomato","turquoise", "violet", "red", "yellowgreen",    "peachpuff", "peru", "pink", "plum2", "purple","wheat", "cornsilk3","cornsilk","coral4","coral",
           "chocolate4","chocolate","black","chartreuse4","chartreuse","burlywood4","burlywood","brown4","blue","bisque4","bisque",
           "azure4","azure","aquamarine4","aquamarine","antiquewhite4","antiquewhite","aliceblue","dodgerblue4","dodgerblue","dimgrey",
           "deepskyblue4", "deepskyblue", "deeppink4", "deeppink","darkviolet", "darkslategray4","darkslategray",
           "darkseagreen4", "darkseagreen", "darksalmon", "darkred", "darkorchid4", "darkorchid", "darkorange4",
           "darkorange","firebrick","darkgreen", "darkgoldenrod4","grey", "darkgoldenrod", "darkcyan")
V(g)$color <- colbar[fun.comps]

set.seed(123)
pdf("C:/graduate/WH赤霉病/pdf picture/16S FZ network module.pdf",height=8,width =9)
plot(g, edge.width=0.5, vertex.frame.color=NA,vertex.label=NA,edge.lty=1,     edge.curved=T,vertex.size=3,margin=c(0, 0,0,0),main="Symptomatic Root",cex.main=0.5) 
dev.off()

da0<- df.tmp4
da0[is.na(da0)] = 0
da <- da0[(da0$kingdom == "Bacteria")& (da0$kingdom.1 == "Bacteria"),]
g <- graph.data.frame(da, directed=FALSE)

M1$BB<-da$color
E(g)$color = as.character(M1$BB)

fun.fc<-cluster_fast_greedy(g)
print(modularity(fun.fc))
modularity(fun.fc,membership(fun.fc))
membership(fun.fc)
print(sizes(fun.fc))###
fun.comps <- membership(fun.fc)
colbar <-c("#0000FF","#FF3399","#FFCC33","#9999FF", "#FF99FF", "#FFFFCC",
           "#33FFFF","#CCFF00","#CCFFCC","#FF3300","#CC99CC","#CC66CC","#FFFFCC","#CC0033","#666666","slateblue1","springgreen","steelblue1","tan1","thistle1","tomato","turquoise", "violet", "red", "yellowgreen",    "peachpuff", "peru", "pink", "plum2", "purple","wheat", "cornsilk3","cornsilk","coral4","coral",
           "chocolate4","chocolate","black","chartreuse4","chartreuse","burlywood4","burlywood","brown4","blue","bisque4","bisque",
           "azure4","azure","aquamarine4","aquamarine","antiquewhite4","antiquewhite","aliceblue","dodgerblue4","dodgerblue","dimgrey",
           "deepskyblue4", "deepskyblue", "deeppink4", "deeppink","darkviolet", "darkslategray4","darkslategray",
           "darkseagreen4", "darkseagreen", "darksalmon", "darkred", "darkorchid4", "darkorchid", "darkorange4",
           "darkorange","firebrick","darkgreen", "darkgoldenrod4","grey", "darkgoldenrod", "darkcyan")
V(g)$color <- colbar[fun.comps]

set.seed(123)
pdf("C:/graduate/WH赤霉病/pdf picture/16S HZ network module.pdf",height=8,width =9)
plot(g, edge.width=0.5, vertex.frame.color=NA,vertex.label=NA,edge.lty=1,     edge.curved=T,vertex.size=3,margin=c(0, 0,0,0),main="Symptomless Root",cex.main=0.5) 
dev.off()

###################Fungal network module###################
setwd("C:/graduate/WH赤霉病/扩增子/WH Fungi")
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
# Module fun-fun###
####################     
FF0<-read.csv("IDotutab.rar.csv", head = T, row.names = 1)
ID.tmp<-FF0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus",	"OTU.ID")]

da0<- df.tmp1
da0[is.na(da0)] = 0
da <- da0[da0$kingdom == "Fungi"&da0$kingdom.1 == "Fungi",]
g <- graph.data.frame(da, directed=FALSE)

M1$FF<-da$color
E(g)$color = as.character(M1$FF)

fun.fc<-cluster_fast_greedy(g)
print(modularity(fun.fc))
modularity(fun.fc,membership(fun.fc))
membership(fun.fc)
print(sizes(fun.fc))###
fun.comps <- membership(fun.fc)
colbar <-c("#0000FF","#FF3399","#FFCC33","#9999FF", "#FF99FF", "#FFFFCC",
           "#33FFFF","#CCFF00","#CCFFCC","#FF3300","#CC99CC","#CC66CC","#FFFFCC","#CC0033","#666666","slateblue1","springgreen","steelblue1","tan1","thistle1","tomato","turquoise", "violet", "red", "yellowgreen",    "peachpuff", "peru", "pink", "plum2", "purple","wheat", "cornsilk3","cornsilk","coral4","coral",
           "chocolate4","chocolate","black","chartreuse4","chartreuse","burlywood4","burlywood","brown4","blue","bisque4","bisque",
           "azure4","azure","aquamarine4","aquamarine","antiquewhite4","antiquewhite","aliceblue","dodgerblue4","dodgerblue","dimgrey",
           "deepskyblue4", "deepskyblue", "deeppink4", "deeppink","darkviolet", "darkslategray4","darkslategray",
           "darkseagreen4", "darkseagreen", "darksalmon", "darkred", "darkorchid4", "darkorchid", "darkorange4",
           "darkorange","firebrick","darkgreen", "darkgoldenrod4","grey", "darkgoldenrod", "darkcyan")
V(g)$color <- colbar[fun.comps]

set.seed(123)
pdf("C:/graduate/WH赤霉病/pdf picture/ITS FR network module.pdf",height=8,width =9)
plot(g, edge.width=0.5, vertex.frame.color=NA,vertex.label=NA,edge.lty=1,     edge.curved=T,vertex.size=3,margin=c(0, 0,0,0),main="Symptomatic Rhizosphere",cex.main=0.5) 
dev.off()


da0<- df.tmp2
da0[is.na(da0)] = 0
da <- da0[da0$kingdom == "Fungi"&da0$kingdom.1 == "Fungi",]
g <- graph.data.frame(da, directed=FALSE)

M1$FF<-da$color
E(g)$color = as.character(M1$FF)

fun.fc<-cluster_fast_greedy(g)
print(modularity(fun.fc))
modularity(fun.fc,membership(fun.fc))
membership(fun.fc)
print(sizes(fun.fc))###
fun.comps <- membership(fun.fc)
colbar <-c("#0000FF","#FF3399","#FFCC33","#9999FF", "#FF99FF", "#FFFFCC",
           "#33FFFF","#CCFF00","#CCFFCC","#FF3300","#CC99CC","#CC66CC","#FFFFCC","#CC0033","#666666","slateblue1","springgreen","steelblue1","tan1","thistle1","tomato","turquoise", "violet", "red", "yellowgreen",    "peachpuff", "peru", "pink", "plum2", "purple","wheat", "cornsilk3","cornsilk","coral4","coral",
           "chocolate4","chocolate","black","chartreuse4","chartreuse","burlywood4","burlywood","brown4","blue","bisque4","bisque",
           "azure4","azure","aquamarine4","aquamarine","antiquewhite4","antiquewhite","aliceblue","dodgerblue4","dodgerblue","dimgrey",
           "deepskyblue4", "deepskyblue", "deeppink4", "deeppink","darkviolet", "darkslategray4","darkslategray",
           "darkseagreen4", "darkseagreen", "darksalmon", "darkred", "darkorchid4", "darkorchid", "darkorange4",
           "darkorange","firebrick","darkgreen", "darkgoldenrod4","grey", "darkgoldenrod", "darkcyan")
V(g)$color <- colbar[fun.comps]

set.seed(123)
pdf("C:/graduate/WH赤霉病/pdf picture/ITS HR network module.pdf",height=8,width =9)
plot(g, edge.width=0.5, vertex.frame.color=NA,vertex.label=NA,edge.lty=1,     edge.curved=T,vertex.size=3,margin=c(0, 0,0,0),main="Symptomless Rhizosphere",cex.main=0.5) 
dev.off()


da0<- df.tmp3
da0[is.na(da0)] = 0
da <- da0[da0$kingdom == "Fungi"&da0$kingdom.1 == "Fungi",]
g <- graph.data.frame(da, directed=FALSE)

M1$FF<-da$color
E(g)$color = as.character(M1$FF)

fun.fc<-cluster_fast_greedy(g)
print(modularity(fun.fc))
modularity(fun.fc,membership(fun.fc))
membership(fun.fc)
print(sizes(fun.fc))###
fun.comps <- membership(fun.fc)
colbar <-c("#0000FF","#FF3399","#FFCC33","#9999FF", "#FF99FF", "#FFFFCC",
           "#33FFFF","#CCFF00","#CCFFCC","#FF3300","#CC99CC","#CC66CC","#FFFFCC","#CC0033","#666666","slateblue1","springgreen","steelblue1","tan1","thistle1","tomato","turquoise", "violet", "red", "yellowgreen",    "peachpuff", "peru", "pink", "plum2", "purple","wheat", "cornsilk3","cornsilk","coral4","coral",
           "chocolate4","chocolate","black","chartreuse4","chartreuse","burlywood4","burlywood","brown4","blue","bisque4","bisque",
           "azure4","azure","aquamarine4","aquamarine","antiquewhite4","antiquewhite","aliceblue","dodgerblue4","dodgerblue","dimgrey",
           "deepskyblue4", "deepskyblue", "deeppink4", "deeppink","darkviolet", "darkslategray4","darkslategray",
           "darkseagreen4", "darkseagreen", "darksalmon", "darkred", "darkorchid4", "darkorchid", "darkorange4",
           "darkorange","firebrick","darkgreen", "darkgoldenrod4","grey", "darkgoldenrod", "darkcyan")
V(g)$color <- colbar[fun.comps]

set.seed(123)
pdf("C:/graduate/WH赤霉病/pdf picture/ITS FZ network module.pdf",height=8,width =9)
plot(g, edge.width=0.5, vertex.frame.color=NA,vertex.label=NA,edge.lty=1,     edge.curved=T,vertex.size=3,margin=c(0, 0,0,0),main="Symptomatic Root",cex.main=0.5) 
dev.off()


da0<- df.tmp4
da0[is.na(da0)] = 0
da <- da0[da0$kingdom == "Fungi"&da0$kingdom.1 == "Fungi",]
g <- graph.data.frame(da, directed=FALSE)

M1$FF<-da$color
E(g)$color = as.character(M1$FF)

fun.fc<-cluster_fast_greedy(g)
print(modularity(fun.fc))
modularity(fun.fc,membership(fun.fc))
membership(fun.fc)
print(sizes(fun.fc))###
fun.comps <- membership(fun.fc)
colbar <-c("#0000FF","#FF3399","#FFCC33","#9999FF", "#FF99FF", "#FFFFCC",
           "#33FFFF","#CCFF00","#CCFFCC","#FF3300","#CC99CC","#CC66CC","#FFFFCC","#CC0033","#666666","slateblue1","springgreen","steelblue1","tan1","thistle1","tomato","turquoise", "violet", "red", "yellowgreen",    "peachpuff", "peru", "pink", "plum2", "purple","wheat", "cornsilk3","cornsilk","coral4","coral",
           "chocolate4","chocolate","black","chartreuse4","chartreuse","burlywood4","burlywood","brown4","blue","bisque4","bisque",
           "azure4","azure","aquamarine4","aquamarine","antiquewhite4","antiquewhite","aliceblue","dodgerblue4","dodgerblue","dimgrey",
           "deepskyblue4", "deepskyblue", "deeppink4", "deeppink","darkviolet", "darkslategray4","darkslategray",
           "darkseagreen4", "darkseagreen", "darksalmon", "darkred", "darkorchid4", "darkorchid", "darkorange4",
           "darkorange","firebrick","darkgreen", "darkgoldenrod4","grey", "darkgoldenrod", "darkcyan")
V(g)$color <- colbar[fun.comps]

set.seed(123)
pdf("C:/graduate/WH赤霉病/pdf picture/ITS HZ network module.pdf",height=8,width =9)
plot(g, edge.width=0.5, vertex.frame.color=NA,vertex.label=NA,edge.lty=1,     edge.curved=T,vertex.size=3,margin=c(0, 0,0,0),main="Symptomless Root",cex.main=0.5) 
dev.off()
