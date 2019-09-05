########HBV progression diversity and venn graph
#######Fig2.A  Venn graph for the four HBV progression group
library(vegan)
library(ggplot2)
install.packages("VennDiagram")
library(VennDiagram)
packageVersion("VennDiagram")
library(gridExtra)
####Set the correct pathway containing the input file
setwd("./dir")
otus_table=read.csv("otus_table.csv",sep=",",header=T,row.names=1)
HBV_env=read.csv("HBV_infection_stage.csv",sep=",",row.names=1,header=T)
all(row.names(HBV_env)==row.names(otus_table))

HBV_otu_env<-cbind(HBV_env,otus_table)
write.csv(HBV_otu_env,file="HBV_otu_env.csv")
#######Subset to each HBV group
Healthy<-HBV_otu_env[HBV_otu_env$Category=="Healthy",]
HBVI<-HBV_otu_env[HBV_otu_env$Category=="HBVI",]
CHB<-HBV_otu_env[HBV_otu_env$Category=="CHB",]
LC<-HBV_otu_env[HBV_otu_env$Category=="LC",]
##remove all columns where all values=0
H_venn<-Healthy[,which(!apply(Healthy,2,FUN = function(x){all(x == 0)}))]
HBVI_venn<-HBVI[,which(!apply(HBVI,2,FUN = function(x){all(x == 0)}))]
CHB_venn<-CHB[,which(!apply(CHB,2,FUN = function(x){all(x == 0)}))]
LC_venn<-LC[,which(!apply(LC,2,FUN = function(x){all(x == 0)}))]

##get the list of vector about the OTU name 
a<-colnames(H_venn[-(1:4)])
b<-colnames(HBVI_venn[-(1:4)])
c<-colnames(CHB_venn[-(1:4)])
d<-colnames(LC_venn[-(1:4)])
list.d=write.csv(c,"c.csv")
###conmbine list.a, list.b, list.c and list.d to obtain list
list=read.csv("list.csv",header=T,sep=",")
venn.diagram(x=list(Healthy=list$Healthy,HBVI=list$HBVI,CHB=list$CHB,LC=list$LC),"venne.tiff",height = 3000, width = 3000, cat.cex=0.9, 
             cex=0.8,margin=c(0.05,0.05,0.05,0.05),col = "transparent",fill = c("#F8766E","#7CAD00", "#00BEC3","#C67CFF"),
             print.mode = "percent",alpha = 0.50)

#########Fig2.B Rarefaction curve observed curve
####please see the details in dirctory Figure2.B

#######Fig2. C Diversity index graph
library(lattice)
library(vegan)
library(tidyr)
load("diversity.total.HBV.rda")
######diversity index among groups
diversity.total.HBV$Category=factor(diversity.total.HBV$Category, levels=c("Healthy", "HBVI", "CHB","LC"), order=TRUE) 
Simpson.plot<-ggplot(diversity.total.HBV,aes(x=Category,y=Simpson,fill=Category))+theme(axis.title.x = element_blank())+
  geom_boxplot(alpha=0.7,width=0.7,outlier.size=0.8)+stat_summary(geom="point",fun.y=mean,pch=8)+
  theme(legend.position="none")+ geom_point()+
  scale_fill_manual(values=c("#F8766E","#7CAD00", "#00BEC3","#C67CFF"))+
  theme(strip.text.x = element_text(size=12,face="bold"),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="black", fill="grey"))+   #modify facet label appearance 
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=10,face="bold"),
        axis.title.y = element_text(face="bold", colour="black", size=15),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=12))+theme(axis.text.x = element_blank())

H.plot<-ggplot(diversity.total.HBV,aes(x=Category,y=Shannon,fill=Category))+theme(axis.title.x = element_blank())+
  geom_boxplot(alpha=0.7,width=0.7,outlier.size=0.8)+stat_summary(geom="point",fun.y=mean,pch=8)+
  theme(legend.position="none")+geom_point()+
  scale_fill_manual(values=c("#F8766E","#7CAD00", "#00BEC3","#C67CFF"))+
  theme(strip.text.x = element_text(size=12,face="bold"),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="black", fill="grey"))+   #modify facet label appearance 
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=10,face="bold"),
        axis.title.y = element_text(face="bold", colour="black", size=15),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=12))

J.plot<-ggplot(diversity.total.HBV,aes(x=Category,y=J,fill=Category))+theme(axis.title.x = element_blank())+
  geom_boxplot(alpha=0.7,width=0.7,outlier.size=0.8)+stat_summary(geom="point",fun.y=mean,pch=8)+
  theme(legend.position="none")+ geom_point()+scale_fill_manual(values=c("#F8766E","#7CAD00", "#00BEC3","#C67CFF"))+
  theme(strip.text.x = element_text(size=12,face="bold"),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="black", fill="grey"))+   #modify facet label appearance 
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=10,face="bold"),
        axis.title.y = element_text(face="bold", colour="black", size=15),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=12))+theme(axis.text.x = element_blank())

chao1.plot<-ggplot(diversity.total.HBV,aes(x=Category,y=Chao1,fill=Category))+theme(axis.title.x = element_blank())+
  geom_boxplot(alpha=0.7,width=0.7,outlier.size=0.8)+stat_summary(geom="point",fun.y=mean,pch=8)+
  theme(legend.position="none")+ geom_point()+scale_fill_manual(values=c("#F8766E","#7CAD00", "#00BEC3","#C67CFF"))+
  theme(strip.text.x = element_text(size=12,face="bold"),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="black", fill="grey"))+   #modify facet label appearance 
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=10,face="bold"),
        axis.title.y = element_text(face="bold", colour="black", size=15),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=12))

######muitiplot function to combine plot together which is list in the following 
multiplot(Simpson.plot,H.plot,J.plot,chao1.plot,cols=2)

#########Group test for Figure 2.C
####Firtst use fligner.test to see whether equal variance
####Then non-parameter method kruskal.test to see whether significant difference among groups
####pairwise.t.test to see how difference among groups
#####Simpson
fligner.test(Simpson~Category,diversity.total.HBV)
kruskal.test(diversity.total.HBV$Simpson~diversity.total.HBV$Category)
ptt.rst <- pairwise.t.test(diversity.total.HBV$Simpson, diversity.total.HBV$Category,p.adjust.method = "BH", pool.sd = FALSE)
grps <- get_groups(ptt.rst, alpha = 0.05, rm.subset = FALSE)
grps$groups
#####Shannon
fligner.test(Shannon~Category,diversity.total.HBV)
kruskal.test(diversity.total.HBV$Shannon~diversity.total.HBV$Category)
ptt.rst <- pairwise.t.test(diversity.total.HBV$Shannon, diversity.total.HBV$Category,p.adjust.method = "BH", pool.sd = FALSE)
grps <- get_groups(ptt.rst, alpha = 0.05, rm.subset = FALSE)
grps$groups
#####J
fligner.test(J~Category,diversity.total.HBV)
kruskal.test(diversity.total.HBV$J~diversity.total.HBV$Category)
ptt.rst <- pairwise.t.test(diversity.total.HBV$J, diversity.total.HBV$Category,p.adjust.method = "BH", pool.sd = FALSE)
grps <- get_groups(ptt.rst, alpha = 0.05, rm.subset = FALSE)
grps$groups
#####Chao1
fligner.test(Chao1~Category,diversity.total.HBV)
kruskal.test(diversity.total.HBV$Chao1~diversity.total.HBV$Category)
ptt.rst <- pairwise.t.test(diversity.total.HBV$Chao1, diversity.total.HBV$Category,p.adjust.method = "BH", pool.sd = FALSE)
grps <- get_groups(ptt.rst, alpha = 0.05, rm.subset = FALSE)
grps$groups

##########Multiplot function
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##########get_groups2
get_groups <- function(paired_t_test_result, alpha = 0.05, rm.subset = FALSE) {
  
  # Get a matrix of the alpha values
  temp <- paired_t_test_result$p.value
  
  # Make a square matrix to populate with the alpha values.
  n <- nrow(temp)
  mat.names <- c(colnames(temp), rownames(temp)[n])
  my.mat <- matrix(data = NA, nrow = n+1, ncol = n+1)
  colnames(my.mat) <- mat.names
  rownames(my.mat) <- mat.names
  
  # Add diagonal.
  for (i in 1:nrow(my.mat)) {
    my.mat[ i, i] <- 0
  }
  
  # Get vector of p.values
  stat <- na.exclude(as.vector(paired_t_test_result$p.value))
  
  # Add other cells to square matrix.
  k=1
  for (j in 1:(nrow(my.mat)-1)) {
    for (i in ((j+1):nrow(my.mat))) {
      my.mat[i,j] <-  my.mat[j,i] <- stat[k]
      k=k+1
    }
  }
  
  # For each column, get list of treatments not significantly different.
  grp <- list()
  trts <- colnames(my.mat)
  for (i in 1:ncol(my.mat)) {
    grp[[i]] <-c(trts[i], names(which(my.mat[ , i] > alpha)))
  }
  
  # Remove groups that are sub-sets of other groups
  k <- 0
  del <- vector()
  for (i in 1:(length(grp)-1)) {
    for ( j in (i+1):length(grp)) {
      if (!rm.subset) {
        if (setequal(grp[[i]], grp[[j]])) {
          k <- k+1
          del[k] <- j
        }
      }
      else {
        if (all(is.element(grp[[i]], grp[[j]]))) {
          k <- k+1
          del[k] <- i
        }
        else if (all(is.element(grp[[j]], grp[[i]]))) {
          k <- k+1
          del[k] <- j
        }
      }
      
    }
  }
  
  del <- unique(del)
  del <- del[order(del, decreasing = TRUE)]
  
  for (i in 1:length(del)) {
    grp[[del[i]]] <- NULL
  }
  
  return(list(groups = grp, p.matrix = my.mat))
  
}
