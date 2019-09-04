############Fig.6
#####Create node and edge input file for gephi 0.9.1
#####Subset the total to health, HBVI, CHB and LC group separately to create four network
#####create node file as Y.cor.mat.test.csv for each group separately
setwd("./dir")
library(psych)
Ycor=read.csv(file.choose(),row.names=1)
attach(Ycor)
Y.cor.test.mat=corr.test(x=Ycor,y=NULL,use = "pairwise",method="pearson",adjust="fdr",alpha=.05,ci=TRUE)
mySTLdata.DF <- as.data.frame(Y.cor.test.mat$r)
write.table(mySTLdata.DF,sep=",",file="Y.cor.mat.test.csv")
#####In the node file, genus name as ID and Phyla name as class 

#####Create edge file as cor_matrix_r.csv for each group
cor_matrix=read.csv(file.choose(),row.names=1)

a<-data.frame(
  Row=rownames(cor_matrix)[row(cor_matrix)[lower.tri(cor_matrix)]],
  Col=colnames(cor_matrix)[col(cor_matrix)[lower.tri(cor_matrix)]],
  Corr=cor_matrix[lower.tri(cor_matrix)])
write.table(a,quote=FALSE,sep=",",file="cor_matrix_r.csv")

#####In the edge file, each genus as Source and the other related with Target and the value following (both positive and negative indicates the 
#####correlation between genus based on their relative abundance within each group, and only
#####absolute value >0.6 were retained.)




