#########Code for Fig.1
#####library necessary R package
library(vegan)
library(ggplot2)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
#####set path to dir containing the input file
setwd("~./dir")
########Fig1.A Ordination plots of different stages with the clinical data
clinical=read.csv("clinical_adjust.csv",sep=",",row.names=1,header=T)
data1=clinical[,-(1:4)][,-(26:27)]
group1=clinical[,1:4]
ord=metaMDS(data1,dist="bray")
group1$Category=factor(group1$Category,levels=c("Healthy","HBVI","CHB","LC"))
scrs <- scores(ord, display = 'sites')
scrs <- cbind(as.data.frame(scrs), Category= group1$Category)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Category, data = scrs, FUN = mean)
segs <- merge(scrs, setNames(cent, c('Category','oNMDS1','oNMDS2')),
              by = 'Category', sort = FALSE)
svg("Fig1.A.svg")
Fig1.A=ggplot(scrs, aes(x = NMDS1, y = NMDS2, colour = Category)) +
  geom_segment(data = segs,
               mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # spiders
  geom_point(data = cent, size = 5) +                         # centroids
  geom_point()+                                              # sample scores
  coord_fixed()+theme(legend.title=element_blank()) +
  theme(axis.title.x = element_text(face="bold", colour="black", size=12),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=12,face="bold"))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=12),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=12,face="bold"))+
  theme(legend.title = element_text(colour="black", size=10, face="bold"))+
  geom_text(x=0.07, y=0.18, label="Stress=0.14",color="black")+
  theme(legend.justification=c(0.05,0.95), legend.position=c(0.05,0.95))
plot(Fig1.A)
dev.off()

#######Fig1.B  db-RDA plots, envfit clinical data to bacterial community
data.adj=read.csv("OTU_adjust.csv",sep=",",row.names=1,header=T)
data.adj.re=decostand(data.adj,"hellinger")
data.dist=vegdist(data.adj.re,method="bray")
clinical.re=decostand(clinical[,-(1:4)],"standardize")

par(mar=c(5,9,1,9))
HBV.dbrda<-capscale(data.dist~WBC+RBC+Hb+PLT+TP+ALB+GLO+ALB.GLO+TBIL+DBIL+IBIL+ALT+AST+AST.ALT+GGTL+ALP+BUN+                     SCR+BUN.SCR+BUA+TG+CHOL+HDL_C+LDL_C+FBG+HbsAg+HBVDNA,clinical.re,dist="bray")
#####First to see which clinical parameter coordinates with bacterial community well and sort them out
fit001<-envfit(HBV.dbrda~WBC+RBC+Hb+PLT+TP+ALB+GLO+ALB.GLO+TBIL+DBIL+IBIL+ALT+AST+AST.ALT+GGTL+ALP+BUN+
                 SCR+BUN.SCR+BUA+TG+CHOL+HDL_C+LDL_C+FBG+HbsAg+HBVDNA,clinical.re)
anova(HBV.dbrda)
#####Second envfit those clinical parameters below to the ordination plots
fit002<-envfit(HBV.dbrda~WBC+RBC+Hb+PLT+TP+ALB+GLO+ALB.GLO+AST.ALT+HDL_C,clinical.re)
plot(HBV.dbrda,type="n",xlab="db-RDA1",ylab="db-RDA2",las=1)
points(HBV.dbrda,dis="si",pch=c(19),col=c("#F8766E","#7CAD00","#00BEC3","#C67CFF"))
plot(fit002,col="darkblue",cex=0.7,font=0.6)
text(x=4.0,y=-1.0,labels="p<0.01")

#######Fig1.C Ordination plots of different stages in Genus levels
######Read otu table files
otus_table=read.csv("otus_table.csv",sep=",",header=T,row.names=1)
otus.relative_abundance=decostand(otus_table,method="hellinger")
write.csv(otus.relative_abundance,"otus.relative_abundance.csv")
######Read group data
HBV_env=read.csv("HBV_infection_stage.csv",sep=",",header=T,row.names=1)
######Check both row names
all(row.names(otus.relative_abundance)==row.names(HBV_env))
######merge all files
total.HBV.Env=cbind(HBV_env,otus.relative_abundance)
#######NMDS analysis for the total communities of six patients
ord.HBV=metaMDS(otus.relative_abundance,distance="bray")
HBV_env$Category = factor(HBV_env$Category, levels=c("Healthy", "HBVI", "CHB","LC")) 
scrs.otu <- scores(ord.HBV, display = 'sites')
scrs.otu <- cbind(as.data.frame(scrs.otu), Category=HBV_env$Category)
cent.otu <- aggregate(cbind(NMDS1, NMDS2) ~ Category, data = scrs.otu, FUN = mean)
segs.otu <- merge(scrs.otu, setNames(cent.otu, c('Category','oNMDS1','oNMDS2')),
              by = 'Category', sort = FALSE)
svg("Fig1.C.svg")
Fig1.C=ggplot(scrs.otu, aes(x = NMDS1, y = NMDS2, colour = Category)) +
  geom_segment(data = segs.otu,
               mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # spiders
  geom_point(data = cent.otu, size = 5) +                         # centroids
  geom_point()+                                              # sample scores
  coord_fixed()+
  theme(axis.title.x = element_text(face="bold", colour="black", size=12),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=12,face="bold"))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=12),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=12,face="bold"))+
  geom_text(x=1.2, y=0.5, label="Stress=0.20",color="black")+
  ggtitle("Bray-Curtis/OTU")+theme(plot.title = element_text(hjust = 0))
plot(Fig1.C)
dev.off()
#######calculate the pairwise Adonis R2 value
pairwise.adonis(otus.relative_abundance,HBV_env$Category,
                sim.method="bray",p.adjust.m="bonferroni")


#######Fig1.D Ordination plots of different stages in Genus levels
genus=read.csv("genus.csv",sep=",",header=T,row.names = 1)
######relative abundance of genus
genus.relative_abundance=decostand(genus,method="hellinger")
#####check both files
all(row.names(HBV_env)==row.names(genus.relative_abundance))
######merge all files
total.HBV.genus=cbind(HBV_env,genus.relative_abundance)
#######NMDS analysis for the total communities of six patients in genus level
ord.genus=metaMDS(genus.relative_abundance,distance="bray")
HBV_env$Category = factor(HBV_env$Category, levels=c("Healthy", "HBVI", "CHB","LC"))  
scrs.genus <- scores(ord.genus, display = 'sites')
scrs.genus <- cbind(as.data.frame(scrs.genus), Category=HBV_env$Category)
cent.genus <- aggregate(cbind(NMDS1, NMDS2) ~ Category, data = scrs.genus, FUN = mean)
segs.genus <- merge(scrs.genus, setNames(cent.genus, c('Category','oNMDS1','oNMDS2')),
                  by = 'Category', sort = FALSE)
svg("Fig1.D.svg")
Fig1.D=ggplot(scrs.genus, aes(x = NMDS1, y = NMDS2, colour = Category)) +
  geom_segment(data = segs.genus,
               mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # spiders
  geom_point(data = cent.genus, size = 5) +                         # centroids
  geom_point()+                                              # sample scores
  coord_fixed()+
  theme(axis.title.x = element_text(face="bold", colour="black", size=12),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=12,face="bold"))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=12),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=12,face="bold"))+
  geom_text(x=0.9, y=0.5, label="Stress=0.22",color="black")+
  ggtitle("Bray-Curtis/Genus")+theme(plot.title = element_text(hjust = 0))
Fig1.D
dev.off()
pairwise.adonis(genus.relative_abundance,HBV_env$Category,
                sim.method="bray",p.adjust.m="bonferroni")


