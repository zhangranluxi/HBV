########Fig.3 
######HBV phyla and class distribution
#####read phyla data
setwd("./dir")

library(vegan)
library(ggplot2)

#######Fig3.A
phyla=read.csv("Phylum_data.csv",sep=",",header=T,row.names=1)
HBV_stage=read.csv("HBV_infection_stage.csv",sep=",",header=T, row.names=1)
all(row.names(phyla)==row.names(HBV_stage))

phyla_rela=decostand(phyla[,-1],method="total")*100

#####remove phyla in low abundance (<1%)
phyla.relative<-phyla_rela[,apply(phyla_rela,2,function(column) all(column>0.00001))]

#####merge phyla and clinical data together
phyla.all=cbind(HBV_stage,phyla.relative)[,-(1:2)][,-2]

phyla.all$Category=factor(phyla.all$Category,levels=c("Healthy","HBVI","CHB","LC"))

######Calculation of data and summary
phyla.all.ad.Ba <- summarySE(phyla.all, measurevar="p_Bacteroidetes",groupvars=c("Category"))

phyla.all.ad.Fi<-summarySE(phyla.all, measurevar="p_Firmicutes",groupvars=c("Category"))

phyla.all.ad.Pr<-summarySE(phyla.all, measurevar="p_Proteobacteria",groupvars=c("Category"))

phyla.all.ad.Ac<-summarySE(phyla.all, measurevar="p_Actinobacteria",groupvars=c("Category"))

p1=ggplot(phyla.all.ad.Ba, aes(x=Category, y=p_Bacteroidetes,group=1)) + 
  geom_errorbar(aes(ymin=p_Bacteroidetes-se, ymax=p_Bacteroidetes+se), width=.1) +
  geom_line()+geom_point()+ylab("p_Bacteroidetes(%)")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())
p1

p2=ggplot(phyla.all.ad.Fi, aes(x=Category, y=p_Firmicutes,group=1)) + 
  geom_errorbar(aes(ymin=p_Firmicutes-se, ymax=p_Firmicutes+se), width=.1) +
  geom_line()+geom_point()+ylab("p_Firmicutes(%)")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())
p2

p3=ggplot(phyla.all.ad.Pr, aes(x=Category, y=p_Proteobacteria,group=1)) + 
  geom_errorbar(aes(ymin=p_Proteobacteria-se, ymax=p_Proteobacteria+se), width=.1) +
  geom_line()+geom_point()+ylab("p_Proteobacteria(%)")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())
p3

p4=ggplot(phyla.all.ad.Ac, aes(x=Category, y=p_Actinobacteria,group=1)) + 
  geom_errorbar(aes(ymin=p_Actinobacteria-se, ymax=p_Actinobacteria+se), width=.1) +
  geom_line()+geom_point()+ylab("p_Actinobacteria(%)")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())
p4

#######Statics analysis of the relative abundance of the taxa among groups
library(vegan)
library("laercio")
bartlett.test(p_Bacteroidetes~Category,phyla.all)
summary(aov((phyla.all$p_Bacteroidetes~phyla.all$Category)))
LDuncan(aov(phyla.all$p_Bacteroidetes~phyla.all$Category),"phyla.all$Category")

bartlett.test(p_Firmicutes~Category,phyla.all)
summary(aov((phyla.all$p_Firmicutes~phyla.all$Category)))
LDuncan(aov(phyla.all$p_Firmicutes~phyla.all$Category),"phyla.all$Category")

fligner.test(phyla.all$p_Proteobacteria~phyla.all$Category,phyla.all)
kruskal.test(phyla.all$p_Proteobacteria~phyla.all$Category,phyla.all)
ptt.rst <- pairwise.t.test(phyla.all$p_Proteobacteria, phyla.all$Category,p.adjust.method = "BH", pool.sd = FALSE)
grps <- get_groups(ptt.rst, alpha = 0.05, rm.subset = FALSE)
grps$groups

fligner.test(phyla.all$p_Actinobacteria~phyla.all$Category,phyla.all)
kruskal.test(phyla.all$p_Actinobacteria~phyla.all$Category,phyla.all)
ptt.rst <- pairwise.t.test(phyla.all$p_Actinobacteria, phyla.all$Category,p.adjust.method = "BH", pool.sd = FALSE)
grps <- get_groups(ptt.rst, alpha = 0.05, rm.subset = FALSE)
grps$groups

######Fig3.B
class_rela_ad=read.csv("class_rela_ad.csv",sep=",",header=T,row.names=1)
all(rownames(class_rela_ad)==rownames(HBV_stage))
class.total<-cbind(HBV_stage,class_rela_ad)
colnames(class.total)

library(tidyr)

class.top12.long<-gather(class.total,class,Relative_abundance,c_Bacteroidia:c_Verrucomicrobiae)
class.top12.long$Category=factor(class.top12.long$Category,levels=c("Healthy","HBVI","CHB","DC"))
class.top12.long$class=factor(class.top12.long$class,levels=c("c_Bacteroidia","c_Clostridia","c_Bacilli","c_Erysipelotrichi",
                                                              "c_Alphaproteobacteria","c_Betaproteobacteria","c_Gammaproteobacteria",
                                                              "c_Deltaproteobacteria","c_Actinobacteria",
                                                              "c_Coriobacteriia","c_Fusobacteriia","c_Verrucomicrobiae"))

class.top12.long$class <- with(class.top12.long,reorder(class,Relative_abundance,FUN = mean, order = FALSE))
class.total.graph<-ggplot(class.top12.long,aes(x=Category,y=Relative_abundance,fill=Category,colour=Category))+geom_jitter(size=1.5)+
  facet_wrap( ~class, ncol=4,scales="free_y")+theme(axis.title.x = element_blank())+ylab("Relative abundance ")+
  theme(axis.text.x= element_text(angle = 45, vjust = 1, hjust=1,size=8,face="bold"))+
  theme(axis.text.y= element_text(size=10,face="bold"))+
  theme(plot.title = element_text(lineheight=.8, face="bold",size=7))+
  theme(strip.text.x = element_text(size=8, angle=0))+ylab("Relative abundance(%)")


##########top12 class distribution
class.all=cbind(HBV_stage,class_rela_ad)[,-(1:2)][,-2]
class.all$Category=factor(class.all$Category,levels=c("Healthy","HBVI","CHB","LC"))
class.all.ad.Ba <- summarySE(class.all, measurevar="c_Bacteroidia",groupvars=c("Category"))
p5=ggplot(class.all.ad.Ba, aes(x=Category, y=c_Bacteroidia,group=1)) + 
  geom_errorbar(aes(ymin=c_Bacteroidia-se, ymax=c_Bacteroidia+se), width=.1) +
  geom_line()+geom_point()+ylab("c_Bacteroidia(%)")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())
p5

bartlett.test(c_Bacteroidia~Category,class.all)
summary(aov((class.all$c_Bacteroidia~class.all$Category)))
LDuncan(aov(class.all$c_Bacteroidia~class.all$Category),"class.all$Category")

#######
class.all.ad.cl <- summarySE(class.all, measurevar="c_Clostridia",groupvars=c("Category"))
p6=ggplot(class.all.ad.cl, aes(x=Category, y=c_Clostridia,group=1)) + 
  geom_errorbar(aes(ymin=c_Clostridia-se, ymax=c_Clostridia+se), width=.1) +
  geom_line()+geom_point()+ylab("c_Clostridia(%)")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())
p6

bartlett.test(c_Clostridia~Category,class.all)
summary(aov((class.all$c_Clostridia~class.all$Category)))
LDuncan(aov(class.all$c_Clostridia~class.all$Category),"class.all$Category")


class.all.ad.Bac <- summarySE(class.all, measurevar="c_Bacilli",groupvars=c("Category"))
p7=ggplot(class.all.ad.Bac, aes(x=Category, y=c_Bacilli,group=1)) + 
  geom_errorbar(aes(ymin=c_Bacilli-se, ymax=c_Bacilli+se), width=.1) +
  geom_line()+geom_point()+ylab("c_Bacilli(%)")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())
p7

fligner.test(class.all$c_Bacilli~class.all$Category,class.all)
kruskal.test(class.all$c_Bacilli~class.all$Category,class.all)
ptt.rst <- pairwise.t.test(class.all$c_Bacilli, class.all$Category,p.adjust.method = "BH", pool.sd = FALSE)
grps <- get_groups(ptt.rst, alpha = 0.05, rm.subset = FALSE)
grps$groups


class.all.ad.Er <- summarySE(class.all, measurevar="c_Erysipelotrichi",groupvars=c("Category"))
p8=ggplot(class.all.ad.Er, aes(x=Category, y=c_Erysipelotrichi,group=1)) + 
  geom_errorbar(aes(ymin=c_Erysipelotrichi-se, ymax=c_Erysipelotrichi+se), width=.1) +
  geom_line()+geom_point()+ylab("c_Erysipelotrichi(%)")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())
p8

bartlett.test(c_Erysipelotrichi~Category,class.all)
summary(aov((class.all$c_Erysipelotrichi~class.all$Category)))
LDuncan(aov(class.all$c_Erysipelotrichi~class.all$Category),"class.all$Category")

class.all.ad.Al <- summarySE(class.all, measurevar="c_Alphaproteobacteria",groupvars=c("Category"))
p9=ggplot(class.all.ad.Al, aes(x=Category, y=c_Alphaproteobacteria,group=1)) + 
  geom_errorbar(aes(ymin=c_Alphaproteobacteria-se, ymax=c_Alphaproteobacteria+se), width=.1) +
  geom_line()+geom_point()+ylab("c_Alphaproteobacteria(%)")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())
p9

fligner.test(class.all$c_Alphaproteobacteria~class.all$Category,class.all)
kruskal.test(class.all$c_Alphaproteobacteria~class.all$Category,class.all)
ptt.rst <- pairwise.t.test(class.all$c_Alphaproteobacteria, class.all$Category,p.adjust.method = "BH", pool.sd = FALSE)
grps <- get_groups(ptt.rst, alpha = 0.05, rm.subset = FALSE)
grps$groups

class.all.ad.Be <- summarySE(class.all, measurevar="c_Betaproteobacteria",groupvars=c("Category"))
p10=ggplot(class.all.ad.Be, aes(x=Category, y=c_Betaproteobacteria,group=1)) + 
  geom_errorbar(aes(ymin=c_Betaproteobacteria-se, ymax=c_Betaproteobacteria+se), width=.1) +
  geom_line()+geom_point()+ylab("c_Betaproteobacteria(%)")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())
p10

bartlett.test(c_Betaproteobacteria~Category,class.all)
summary(aov((class.all$c_Betaproteobacteria~class.all$Category)))
LDuncan(aov(class.all$c_Betaproteobacteria~class.all$Category),"class.all$Category")


class.all.ad.Ga <- summarySE(class.all, measurevar="c_Gammaproteobacteria",groupvars=c("Category"))
p11=ggplot(class.all.ad.Ga, aes(x=Category, y=c_Gammaproteobacteria,group=1)) + 
  geom_errorbar(aes(ymin=c_Gammaproteobacteria-se, ymax=c_Gammaproteobacteria+se), width=.1) +
  geom_line()+geom_point()+ylab("c_Gammaproteobacteria(%)")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())
p11

fligner.test(class.all$c_Gammaproteobacteria~class.all$Category,class.all)
kruskal.test(class.all$c_Gammaproteobacteria~class.all$Category,class.all)
ptt.rst <- pairwise.t.test(class.all$c_Gammaproteobacteria, class.all$Category,p.adjust.method = "BH", pool.sd = FALSE)
grps <- get_groups(ptt.rst, alpha = 0.05, rm.subset = FALSE)
grps$groups

class.all.ad.De <- summarySE(class.all, measurevar="c_Deltaproteobacteria",groupvars=c("Category"))
p12=ggplot(class.all.ad.De, aes(x=Category, y=c_Deltaproteobacteria,group=1)) + 
  geom_errorbar(aes(ymin=c_Deltaproteobacteria-se, ymax=c_Deltaproteobacteria+se), width=.1) +
  geom_line()+geom_point()+ylab("c_Deltaproteobacteria(%)")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())
p12

bartlett.test(c_Deltaproteobacteria~Category,class.all)
summary(aov((class.all$c_Deltaproteobacteria~class.all$Category)))
LDuncan(aov(class.all$c_Deltaproteobacteria~class.all$Category),"class.all$Category")

class.all.ad.Ac <- summarySE(class.all, measurevar="c_Actinobacteria",groupvars=c("Category"))
p13=ggplot(class.all.ad.Ac, aes(x=Category, y=c_Actinobacteria,group=1)) + 
  geom_errorbar(aes(ymin=c_Actinobacteria-se, ymax=c_Actinobacteria+se), width=.1) +
  geom_line()+geom_point()+ylab("c_Actinobacteria(%)")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(axis.text.x= element_text(angle = 0, vjust = 0.5, hjust=0.5,size=7,face="bold"))+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())+
  theme(axis.text.x= element_text(angle = 0, vjust = 0.5, hjust=0.5,size=7,face="bold"))
p13

fligner.test(class.all$c_Actinobacteria~class.all$Category,class.all)
kruskal.test(class.all$c_Actinobacteria~class.all$Category,class.all)
ptt.rst <- pairwise.t.test(class.all$c_Actinobacteria, class.all$Category,p.adjust.method = "BH", pool.sd = FALSE)
grps <- get_groups(ptt.rst, alpha = 0.05, rm.subset = FALSE)
grps$groups

class.all.ad.Co <- summarySE(class.all, measurevar="c_Coriobacteriia",groupvars=c("Category"))
p14=ggplot(class.all.ad.Co, aes(x=Category, y=c_Coriobacteriia,group=1)) + 
  geom_errorbar(aes(ymin=c_Coriobacteriia-se, ymax=c_Coriobacteriia+se), width=.1) +
  geom_line()+geom_point()+ylab("c_Coriobacteriia(%)")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())+
  theme(axis.text.x= element_text(angle = 0, vjust = 0.5, hjust=0.5,size=7,face="bold"))
p14

bartlett.test(c_Coriobacteriia~Category,class.all)
summary(aov((class.all$c_Coriobacteriia~class.all$Category)))
LDuncan(aov(class.all$c_Coriobacteriia~class.all$Category),"class.all$Category")

class.all.ad.Fu <- summarySE(class.all, measurevar="c_Fusobacteriia",groupvars=c("Category"))
p15=ggplot(class.all.ad.Fu, aes(x=Category, y=c_Fusobacteriia,group=1)) + 
  geom_errorbar(aes(ymin=c_Fusobacteriia-se, ymax=c_Fusobacteriia+se), width=.1) +
  geom_line()+geom_point()+ylab("c_Fusobacteriia(%)")+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())+
  theme(axis.text.x= element_text(angle = 0, vjust = 0.5, hjust=0.5,size=7,face="bold"))
p15

bartlett.test(c_Fusobacteriia~Category,class.all)
summary(aov((class.all$c_Fusobacteriia~class.all$Category)))
LDuncan(aov(class.all$c_Fusobacteriia~class.all$Category),"class.all$Category")

class.all.ad.Ve <- summarySE(class.all, measurevar="c_Verrucomicrobiae",groupvars=c("Category"))
p16=ggplot(class.all.ad.Ve, aes(x=Category, y=c_Verrucomicrobiae,group=1)) + 
  geom_errorbar(aes(ymin=c_Verrucomicrobiae-se, ymax=c_Verrucomicrobiae+se), width=.1) +
  geom_line()+geom_point()+ylab("c_Verrucomicrobiae(%)")+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(), axis.title.x = element_blank())+
  theme(axis.text.x= element_text(angle = 0, vjust = 0.5, hjust=0.5,size=7,face="bold"))+
  theme(panel.grid.major = element_blank())+theme(panel.grid.minor = element_blank())
p16

bartlett.test(c_Verrucomicrobiae~Category,class.all)
summary(aov((class.all$c_Verrucomicrobiae~class.all$Category)))
LDuncan(aov(class.all$c_Verrucomicrobiae~class.all$Category),"class.all$Category")

multiplot(p1,p5,p8,p16,p2,p6,p9,p15,p3,p10,p11,p13,p4,p12,p7,p14,cols=4)

#########Fig.3 C heatmap of top genus distribution
###genus
library(pheatmap)
annotate1=read.csv("./annotate1.csv", header=T, row.names=1)
ann_color = list(Category=c(Healthy="#F8766E",HBVI="#7CAD00", CHB="#00BEC3", LC="#C67CFF"))
annotate1$Category = factor(annotate1$Category,levels=c("Healthy", "HBVI", "CHB","LC"))
genus=read.csv("top40genus.csv",sep=",",header=T,row.names=1)
genus.re=t(genus[,-1])
genus.re=t(decostand(genus[,-1],method="hellinger"))

pheatmap(genus.re, scale="row", 
         cluster_cols = F,cluster_rows = F,
         annotation_col =annotate1 ,legend_breaks = c(-3,0,3),legend_labels = c(">=-3","0","<=3"),
         annotation_legend = T,
         annotation_colors = ann_color,
         main="Heatmap of top genus taxa",
         show_colnames = F,
         #number_color=black,
         cellwidth = 4,
         cellheight = 6,
         border_color ="white",
         fontsize_col=6,
         fontsize_row= 6,
         color = colorRampPalette(c("blue", "white", "red"))(100))

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
