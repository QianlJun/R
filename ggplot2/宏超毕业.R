setwd("E:/D/RLPC/R/workdata")
p<-read.table("kegg1.csv",header = T,sep=",")
library(ggplot2)
a<-ggplot(p,aes(x=gene_ratio,y=Description,size=Gene_DE,fill=p_value))
a+geom_point(shape=21,color="black")
b<-qplot(gene.ratio,Description,colour=p_value,size=Gene_DE,data=p)
b+theme(axis.text=element_text(size=8,family="Arial"),axis.title=element_text(size=8,family="Arial"),legend.text=element_text(size=5,family="Arial"),legend.title=element_text(size=5,family="Arial")) + scale_colour_gradientn( limits=c(0,0.05),colours=c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#0072B2","#56B4E9","#CC79A7") ) + xlim(0,3) + labs(y="")+scale_size("count",range=c(1,2))
library('showtext')
showtext_auto(enable = TRUE)
font_add("myfont","Arial.ttf")
warnings()

b+theme(axis.text=element_text(size=8),axis.title=element_text(size=8),
        legend.text=element_text(size=5),legend.title=element_text(size=5)) 
      + scale_colour_gradientn( limits=c(0,0.05),
                                colours=c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#0072B2","#56B4E9","#CC79A7") )
      + xlim(0,3) + labs(y="")+scale_size("count",range=c(1,2))
windowsFonts(myFont = windowsFont("Arial")) 
c<-b+theme(axis.text=element_text(size=8,family="myFont"),axis.title=element_text(size=8,family="myFont"),legend.text=element_text(size=8,family="myFont"),legend.title=element_text(size=8,family="myFont")) + scale_colour_gradientn( limits=c(0,0.05),colours=c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#0072B2","#56B4E9","#CC79A7") ) + xlim(0,3) + labs(y="")+scale_size("count",range=c(4,8))
svg(args[2],w=3,h=3)
Cairo(800,400,file=paste(c, ".svg", sep=""),type="svg",bg="transparent",pointsize=8, units="px",dpi=600)
svg(file = "Railways.svg", width = 3, height = 3)
ggsave(c, file="ratings.pdf", width=6, height=6)
c
p<-read.table("kegg1.csv",header = T,sep=",")
b<-qplot(gene.ratio,Description,colour=p_value,size=Gene_DE,data=p)
c<-b+theme(axis.text=element_text(size=10,family="myFont"),axis.title=element_text(size=10,family="myFont"),legend.text=element_text(size=10,family="myFont"),legend.title=element_text(size=10,family="myFont")) + scale_colour_gradientn( limits=c(0,0.05),colours=c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#0072B2","#56B4E9","#CC79A7") ) + xlim(0,3) + labs(y="")+scale_size("count",range=c(5,10))+theme(legend.key.width=unit(0.5,'cm'))
d<-c+xlab("gene ratio")
ggsave(d, file="ratings.pdf", width=5, height=5)


library(svglite)
ggsave(d, file="ratings.svg", width=5, height=5)

library(ggplot2)
library(svglite)
font_add("Arial","Arial.ttf")
p<-read.table("kegg1.csv",header = T,sep=",")
b<-qplot(gene.ratio,Description,colour=p_value,size=Gene_DE,data=p)
c<-b+theme(axis.text=element_text(size=10,family="Arial",color="Black"),axis.title=element_text(size=10,family="Arial"),legend.text=element_text(size=10,family="Arial"),legend.title=element_text(size=10,family="Arial")) + scale_colour_gradientn( limits=c(0,0.05),colours=c("#D55E00", "#E69F00", "#F0E442", "#009E73", "#0072B2","#56B4E9","#CC79A7") ) + xlim(0,3) + labs(y="")+scale_size("count",range=c(5,10))+theme(legend.key.width=unit(0.5,'cm'))
d<-c+xlab("gene ratio")
ggsave(d, file="ratings1.svg", width=5, height=5)


