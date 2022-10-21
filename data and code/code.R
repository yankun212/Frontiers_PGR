

###Figure.1E~H The difference of bacterial β-diversity (E, F, G, H) and Table.S1
rm(list=ls())
library(ggthemes)
library(ggplot2)
library(factoextra)
library(ggsci)
library(vegan)
genus_table<-read.table("C:/Users/maple/Desktop/file/genus_table.txt",header=T)

#######
###Fig.1E
col<-c("#3c5488","#e64b35","#17b3b5","#5dc2b2")
pca <- prcomp(genus_table[-1])
p1<-fviz_pca_ind(pca,geom.ind="point",pointsize=3,
             col.ind=genus_table$group,
             palette="joo",
             addEllipses = T,mean.point=F,
             ellipse.type="confidence",
             ellipse.level=0.95,
             select.var = list(contrib = 10),
             )+theme_few()+
  theme(text = element_text(size = 15,face="bold",family = "serif"))+
  scale_fill_manual(values = col)+
  scale_color_manual(values = col)

###Adonis
distance <- vegdist(genus_table[-1], method = 'bray')
ub_dis_table <- as.dist(distance, diag = FALSE, upper = FALSE)
adonis(ub_dis_table~group, data=genus_table, permutations = 9999) 

###Anosim 
summary(anosim(distance,genus_table$group,permutations = 999))

####MRPP
mrpp(genus_table[-1],genus_table$group)

#####Fig.1F
pca <- prcomp(genus_table[c(1:12,19:21),][-1])
p2<-fviz_pca_ind(pca,geom.ind="point",pointsize=3,
                 col.ind=genus_table[c(1:12,19:21),]$group,
                 palette="joo",
                 addEllipses = T,mean.point=F,
                 ellipse.type="confidence",
                 ellipse.level=0.95,
                 select.var = list(contrib = 10),
                 )+theme_few()+theme(text = element_text(size = 15,face="bold",family = "serif"))

###Adonis
distance <- vegdist(genus_table[c(1:12,19:21),][-1], method = 'bray')
ub_dis_table <- as.dist(distance, diag = FALSE, upper = FALSE)
adonis(ub_dis_table~group, data=genus_table[c(1:12,19:21),], permutations = 9999) 

###Anosim 
summary(anosim(distance,genus_table[c(1:12,19:21),]$group,permutations = 999))

####MRPP
mrpp(genus_table[c(1:12,19:21),][-1],genus_table[c(1:12,19:21),]$group)

##Fig.1G
pca <- prcomp(genus_table[c(7:21),][-1])
p3<-fviz_pca_ind(pca,geom.ind="point",pointsize=3,
                 col.ind=genus_table[c(7:21),]$group,
                 palette="joo",
                 addEllipses = T,mean.point=F,
                 ellipse.type="confidence",
                 ellipse.level=0.95,
                 select.var = list(contrib = 10),
)+theme_few()+theme(text = element_text(size = 15,face="bold",family = "serif"))

###Adonis
distance <- vegdist(genus_table[c(7:21),][-1], method = 'bray')
ub_dis_table <- as.dist(distance, diag = FALSE, upper = FALSE)
adonis(ub_dis_table~group, data=genus_table[c(7:21),], permutations = 9999) 

###Anosim 
summary(anosim(distance,genus_table[c(7:21),]$group,permutations = 999))

####MRPP
mrpp(genus_table[c(7:21),][-1],genus_table[c(7:21),]$group)

##Fig.1H
pca <- prcomp(genus_table[c(1:6,13:21),][-1])
p4<-fviz_pca_ind(pca,geom.ind="point",pointsize=3,
                 col.ind=genus_table[c(1:6,13:21),]$group,
                 palette="joo",
                 addEllipses = T,mean.point=F,
                 ellipse.type="confidence",
                 ellipse.level=0.95,
                 select.var = list(contrib = 10),
)+theme_few()+theme(text = element_text(size = 15,face="bold",family = "serif"))

###Adonis
distance <- vegdist(genus_table[c(1:6,13:21),][-1], method = 'bray')
ub_dis_table <- as.dist(distance, diag = FALSE, upper = FALSE)
adonis(ub_dis_table~group, data=genus_table[c(1:6,13:21),], permutations = 9999) 

###Anosim 
summary(anosim(distance,genus_table[c(1:6,13:21),]$group,permutations = 999))

####MRPP
mrpp(genus_table[c(1:6,13:21),][-1],genus_table[c(1:6,13:21),]$group)


#####Fig.2A

library(reshape)
library(ggplot2)
library(ggthemes)
phylum_table<-read.table("C:/Users/maple/Desktop/file/Top10_phylum.txt",header=T)

###color
c<-c("#f94040","#90bff9","#41f0ae","#099963","#f7abe8","#ffff00","#de8bf9","#caffca","#b856d7","#ffe0c0","#a0a0a4")

phylum_table_m<-melt(phylum_table)
#### Order by abundance
phylum_table_m$phylum<-factor(phylum_table_m$phylum,levels =as.character(phylum_table$phylum))

ggplot(phylum_table_m, aes(x=variable,y=value, fill=phylum)) +
  geom_bar(stat = "identity",width = 0.6) +
  scale_fill_manual(values = c)+
  theme_few()+
  theme(text = element_text(size = 15,face="bold",family = "serif"))


#####Fig.3 This part of the analysis was performed using the website http://huttenhower.sph.harvard.edu/galaxy/.


########
#####Fig.4 The co-occurrence networks (A) and hub nodes (B)

library(WGCNA)
library(igraph)
library(ggplot2)

table<-read.table("C:/Users/maple/Desktop/file/Species_table.txt",header=T)
###Correlation matrix calculation
#####It is necessary to remove the taxa with more 0 values in the group before analysis
occor<-corAndPvalue(t(table),method="pearson",use="p")
occor.r = occor$cor 
occor.p = occor$p 
occor.r[occor.p>0.01|abs(occor.r)<0.8] = 0
diag(occor.r)<-0
occor.r[is.na(occor.r)]<-0
####edge
sum(abs(occor.r)>0)/2 
####node
sum(colSums(abs(occor.r))>0)
#####Export the data and plot using Gephi
write.csv(occor.r,file="network.csv")

######Fig.4B hub nodes

occor.r[abs(occor.r)>0]=1
adjacency_unweight<-occor.r

## Obtain an undirected network with no weights
igraph <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)
igraph 

##computational node
V(igraph)$degree <- degree(igraph)

set.seed(123)
V(igraph)$modularity <- membership(cluster_fast_greedy(igraph))

###integrated data
nodes_list <- data.frame(
  nodes_id = V(igraph)$name, 
  degree = V(igraph)$degree, 
  modularity = V(igraph)$modularity
)
head(nodes_list) 

####Calculate within-module connectivity (Zi) and among-module connectivity (Pi)
source('hub_node.r')

row.names(nodes_list)<-nodes_list$nodes_id
nodes_list<-nodes_list[-1]

zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)

zi_pi <- na.omit(zi_pi)   # remove NA
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type,shape=type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)




#####Fig.5
library(psych)
library(corrplot)
library(reshape)
library(ggplot2)
library(ggsci)
library(ggthemes)

####Importing environmental data
env<-read.table("C:/Users/maple/Desktop/file/metadata.txt",header=T,sep = "\t")[c(3:14)]
#####genus_table   Only genera with abundance greater than 0.01 were included.
r<-corr.test(genus_table,env)
####r<-corr.test(genus_table,env,adjust = "none")
p1<-corrplot(r$r, insig = "blank",method = "square",pch.col = "black",tl.cex=1,tl.col = "black")

#######Plotting with Means
#######genus_m is the file after refactoring with the melt function
p2<-ggplot(genus_m,aes(x=variable,y=genus,color=variable))+
  geom_point(aes(size=value))+
  scale_size_area(max_size =12)+
  theme_few()+
  theme(text = element_text(size = 15,face="bold",family = "serif"))+
  scale_color_npg()


############Fig.6
library(TITAN2)

glades.titan<- titan(genus_table[-1],env$NO3.N,
                       minSplt = 3, numPerm = 250, boot = TRUE, nBoot = 500, imax = FALSE,
                       ivTot = FALSE, pur.cut = 0.85, rel.cut = 0.85, ncpus = 1, memory = FALSE
)

plot_taxa_ridges(glades.titan2, ytxt.sz=8)


############Fig.7 Community assembly analysis  NTI and BNTI

library(picante)
##### Import ASV tables and trees
comm<-read.table("C:/Users/maple/Desktop/file/ASV_table.txt",header=T) 
phy<-read.tree("C:/Users/maple/Desktop/file/rooted_tree.tre")
comm<-t(comm)
prune_tree<-prune.sample(comm,phy)
phydist<-cophenetic(prune_tree)
####Calculate the NTI
mntd<-ses.mntd(comm,phydist,null.model="taxa.labels",abundance.weighted=T, runs=999)
####The observed βMNTD
comdist.resultS12<-comdistnt(comm,phydist)
###nullmodel
f<-function(){
  g<-randomizeMatrix(comm,null.model=c("frequency","richness","independentswap","trialswap"),iterations=1000)
  fc<-comdist.result<-comdistnt(g,phydist)
}
mt<-mean(replicate(999, f()))
mt.sd<-sd(replicate(999, f()))
BNTI=(comdist.resultS12-mt)/mt.sd

