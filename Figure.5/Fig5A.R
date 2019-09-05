#########Heatmap analysis between indicator species and clinical parameters
########Fig5.A
setwd("~/Dir")
######Correlation between species and clinical parameters
clinical_adjust=read.csv("clinical_adjust.csv",sep=",",header=T,row.names=1)
indicator_species1=read.csv("indicator_species1.csv",sep=",",header=T,row.names=1)
all(row.names(clinical_adjust)==row.names(indicator_species1))

######Significance analysis
sig3=data.frame(matrix(nrow=100,ncol=100))
sig4=data.frame(matrix(nrow=100,ncol=100))

######Correlation 
for(c in 1:10) { for (e in 1:49){
  sig3[c,e]=round(cor(clinical_adjust[,c],indicator_species1[,e]),2)}
}

colnames(sig3)=colnames(indicator_species1)

write.csv(sig3,file="sig3.csv")

######significance
for(c in 1:10) { for (e in 1:49){
  sig4[c,e]=round(cor.test(clinical_adjust[,c],indicator_species1[,e])[[3]],2)}
}
colnames(sig4)=colnames(indicator_species1)
write.csv(sig4,file="sig4.csv")

#####correlation pheatmap
library(pheatmap)
library(grid)
sig3=t(read.csv("sig3.csv",sep=",",header=T,row.names=1))
sig4=t(read.csv("sig4.csv",sep=",",header=T,row.names=1))
sig4[as.matrix(sig4<=0.01)]<-"**"
sig4[as.matrix(sig4)>0.01&as.matrix(sig4)<=0.05]<-"*"                           
sig4[as.matrix(sig4)>0.05]<-""

tiff("correlation_clinics.tiff", width=2400, height=1200,res=300)
pheatmap(t(sig3), scale="row",display_numbers=t(sig4), 
         cluster_cols = F,cluster_rows = F,
         annotation_legend = T,
         main="Correlation between indicator species and clinic",
         show_colnames = T,
         #number_color=black,
         cellwidth = 8,
         cellheight = 8,
         border_color ="white",
         fontsize_col=4,
         fontsize_row=4,
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

########The code of Fig5.B is in Fig4 Directory

########Fig5.C is done by Graphpad Prism 5, the input file is indicator.csv file



