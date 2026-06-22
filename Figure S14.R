##################Cross-kingdom positive network################
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
# significant positive correlations #
#####################################
BF0<-read.csv("Bac-Fun.ID.csv", head = T, row.names = 1)
ID.tmp<-BF0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus","OTU.ID")]
r.cutoff = 0.6
p.cutoff = 0.05

spman.r0 <- spman.d1
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[upper.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[upper.tri(Cor)]], Cor=Cor[upper.tri(Cor)])

P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[upper.tri(P0)]], 
                 col=colnames(P0)[col(P0)[upper.tri(P0)]], p=P0[upper.tri(P0)])

df <- data.frame(Cor.df,  P.df, Compartment = "Rhizosphere", State = "Symptomatic")
da.tmp<-df.sig<- df[df$Cor > r.cutoff & df$p < p.cutoff,] 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp1<-data.frame(da.tmp, M1, M2)

spman.r0 <- spman.d2
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[upper.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[upper.tri(Cor)]], Cor=Cor[upper.tri(Cor)])

P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[upper.tri(P0)]], 
                 col=colnames(P0)[col(P0)[upper.tri(P0)]], p=P0[upper.tri(P0)])

df <- data.frame(Cor.df,  P.df, Compartment = "Rhizosphere", State = "Symptomless")
da.tmp<-df.sig<- df[ df$Cor > r.cutoff & df$p < p.cutoff,] 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp2<-data.frame(da.tmp, M1, M2)


spman.r0 <- spman.d3
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[upper.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[upper.tri(Cor)]], Cor=Cor[upper.tri(Cor)])
P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[upper.tri(P0)]], 
                 col=colnames(P0)[col(P0)[upper.tri(P0)]], p=P0[upper.tri(P0)])

df <- data.frame(Cor.df,  P.df, Compartment = "Root", State = "Symptomatic")
da.tmp<-df.sig<- df[ df$Cor > r.cutoff & df$p < p.cutoff,] 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp3<-data.frame(da.tmp, M1, M2)


spman.r0 <- spman.d4
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[upper.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[upper.tri(Cor)]], Cor=Cor[upper.tri(Cor)])

P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[upper.tri(P0)]], 
                 col=colnames(P0)[col(P0)[upper.tri(P0)]], p=P0[upper.tri(P0)])
df <- data.frame(Cor.df,  P.df, Compartment = "Root", State = "Symptomless")
da.tmp<-df.sig<- df[ df$Cor > r.cutoff & df$p < p.cutoff,] 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp4<-data.frame(da.tmp, M1, M2) 


####################
# Module bac-fun###
####################   
BF0<-read.csv("Diseased Rhizosphere Bac-Fun.ID.csv", head = T, row.names = 1)
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

M1$BF<-"red"
M1$BF[M1$kingdom=="Fungi" & M2$kingdom=="Fungi" ]<-"#70c17f"
M1$BF[(M1$kingdom=="Bacteria" )& (M2$kingdom=="Bacteria" )]<-"#7ca9cc"

E(g)$color = as.character(M1$BF)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/pdf picture/cross-kingdom FR positive network.pdf",height=8,width =12)
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
       legend=c("BB","FF","BF"),
       lty=1.5,col=c("#7ca9cc","#70c17f","red"),
       text.col = c("#7ca9cc","#70c17f","red"),
       cex=1, bty="n")
dev.off()

BF0<-read.csv("Healthy Rhizosphere Bac-Fun.ID.csv", head = T, row.names = 1)
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

M1$BF<-"red"
M1$BF[M1$kingdom=="Fungi" & M2$kingdom=="Fungi" ]<-"#70c17f"
M1$BF[(M1$kingdom=="Bacteria" )& (M2$kingdom=="Bacteria" )]<-"#7ca9cc"

E(g)$color = as.character(M1$BF)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/pdf picture/cross-kingdom HR positive network.pdf",height=8,width =12)
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
       legend=c("BB","FF","BF"),
       lty=1.5,col=c("#7ca9cc","#70c17f","red"),
       text.col = c("#7ca9cc","#70c17f","red"),
       cex=1, bty="n")
dev.off()

BF0<-read.csv("Diseased Root Bac-Fun.ID.csv", head = T, row.names = 1)
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

M1$BF<-"red"
M1$BF[M1$kingdom=="Fungi" & M2$kingdom=="Fungi" ]<-"#70c17f"
M1$BF[(M1$kingdom=="Bacteria" )& (M2$kingdom=="Bacteria" )]<-"#7ca9cc"

E(g)$color = as.character(M1$BF)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/pdf picture/cross-kingdom FZ positive network.pdf",height=8,width =12)
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
       legend=c("BB","FF","BF"),
       lty=1.5,col=c("#7ca9cc","#70c17f","red"),
       text.col = c("#7ca9cc","#70c17f","red"),
       cex=1, bty="n")
dev.off()


BF0<-read.csv("Healthy Root Bac-Fun.ID.csv", head = T, row.names = 1)
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

M1$BF<-"red"
M1$BF[M1$kingdom=="Fungi" & M2$kingdom=="Fungi" ]<-"#70c17f"
M1$BF[(M1$kingdom=="Bacteria" )& (M2$kingdom=="Bacteria" )]<-"#7ca9cc"

E(g)$color = as.character(M1$BF)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/pdf picture/cross-kingdom HZ positive network.pdf",height=8,width =12)
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
       legend=c("BB","FF","BF"),
       lty=1.5,col=c("#7ca9cc","#70c17f","red"),
       text.col = c("#7ca9cc","#70c17f","red"),
       cex=1, bty="n")
dev.off()

######################Cross-kingdom negative network####################
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
# significant negative correlations #
#####################################
BF0<-read.csv("Bac-Fun.ID.csv", head = T, row.names = 1)
ID.tmp<-BF0[,c("kingdom",	"phylum",	"class",	"order",	"family",	"genus","OTU.ID")]
r.cutoff = 0.6
p.cutoff = 0.05

spman.r0 <- spman.d1
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[upper.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[upper.tri(Cor)]], Cor=Cor[upper.tri(Cor)])

P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[upper.tri(P0)]], 
                 col=colnames(P0)[col(P0)[upper.tri(P0)]], p=P0[upper.tri(P0)])

df <- data.frame(Cor.df,  P.df, Compartment = "Rhizosphere", State = "Symptomatic")
da.tmp<-df.sig<- df[df$Cor < -r.cutoff & df$p < p.cutoff,] 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp1<-data.frame(da.tmp, M1, M2)


spman.r0 <- spman.d2
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[upper.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[upper.tri(Cor)]], Cor=Cor[upper.tri(Cor)])

P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[upper.tri(P0)]], 
                 col=colnames(P0)[col(P0)[upper.tri(P0)]], p=P0[upper.tri(P0)])

df <- data.frame(Cor.df,  P.df, Compartment = "Rhizosphere", State = "Symptomless")
da.tmp<-df.sig<- df[ df$Cor < -r.cutoff & df$p < p.cutoff,] 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp2<-data.frame(da.tmp, M1, M2)


spman.r0 <- spman.d3
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[upper.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[upper.tri(Cor)]], Cor=Cor[upper.tri(Cor)])
P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[upper.tri(P0)]], 
                 col=colnames(P0)[col(P0)[upper.tri(P0)]], p=P0[upper.tri(P0)])

df <- data.frame(Cor.df,  P.df, Compartment = "Root", State = "Symptomatic")
da.tmp<-df.sig<- df[ df$Cor < -r.cutoff & df$p < p.cutoff,] 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp3<-data.frame(da.tmp, M1, M2)


spman.r0 <- spman.d4
Cor<-as.matrix(spman.r0$r)
Cor.df<-data.frame(row=rownames(Cor)[row(Cor)[upper.tri(Cor)]], 
                   col=colnames(Cor)[col(Cor)[upper.tri(Cor)]], Cor=Cor[upper.tri(Cor)])

P0<-as.matrix(spman.r0$p)
P.df<-data.frame(row=rownames(P0)[row(P0)[upper.tri(P0)]], 
                 col=colnames(P0)[col(P0)[upper.tri(P0)]], p=P0[upper.tri(P0)])
df <- data.frame(Cor.df,  P.df, Compartment = "Root", State = "Symptomless")
da.tmp<-df.sig<- df[ df$Cor < -r.cutoff & df$p < p.cutoff,] 
V1<-data.frame("v1"=da.tmp$row); V2<-data.frame("v2"=da.tmp$col)
IDsub1<-ID.tmp[ID.tmp$OTU.ID %in% V1$v1, ]; IDsub2<-ID.tmp[ID.tmp$OTU.ID %in% V2$v2, ]
V1$id  <- 1:nrow(V1); V2$id  <- 1:nrow(V2)
M1<-merge(V1, IDsub1, by.x = "v1", by.y = "OTU.ID", all.x= T); M1<-M1[order(M1$id), ]
M2<-merge(V2, IDsub2, by.x = "v2", by.y = "OTU.ID", all.x = T); M2<-M2[order(M2$id), ]
df.tmp4<-data.frame(da.tmp, M1, M2) 


####################
# Module bac-fun###
####################   
BF0<-read.csv("Diseased Rhizosphere N Bac-Fun.ID.csv", head = T, row.names = 1)
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

M1$BF<-"red"
M1$BF[M1$kingdom=="Fungi" & M2$kingdom=="Fungi" ]<-"#70c17f"
M1$BF[(M1$kingdom=="Bacteria" )& (M2$kingdom=="Bacteria" )]<-"#7ca9cc"

E(g)$color = as.character(M1$BF)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/pdf picture/cross-kingdom FR negative network.pdf",height=8,width =12)
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
       legend=c("BB","FF","BF"),
       lty=1.5,col=c("#7ca9cc","#70c17f","red"),
       text.col = c("#7ca9cc","#70c17f","red"),
       cex=1, bty="n")
dev.off()

BF0<-read.csv("Healthy Rhizosphere N Bac-Fun.ID.csv", head = T, row.names = 1)
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

M1$BF<-"red"
M1$BF[M1$kingdom=="Fungi" & M2$kingdom=="Fungi" ]<-"#70c17f"
M1$BF[(M1$kingdom=="Bacteria" )& (M2$kingdom=="Bacteria" )]<-"#7ca9cc"

E(g)$color = as.character(M1$BF)

set.seed(123)
par(mfrow=c(1,1),mar=c(0, 0, 0, 0))
pdf("C:/graduate/WH赤霉病/pdf picture/cross-kingdom HR negative network.pdf",height=8,width =12)
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
       legend=c("BB","FF","BF"),
       lty=1.5,col=c("#7ca9cc","#70c17f","red"),
       text.col = c("#7ca9cc","#70c17f","red"),
       cex=1, bty="n")
dev.off()
