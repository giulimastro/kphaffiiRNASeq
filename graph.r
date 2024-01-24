library("ggplot2")
library("DESeq2")
library("reshape2")
library("stringr")

vsd <- assay(normTransform(dds))

plasmids<-data.frame(tail(vsd,9))
plasmids$plasmid=rownames(plasmids)

plasmids<-melt(plasmids)

variables=str_split_fixed(plasmids$variable,"\\.",3)
colnames(variables)<-c("VLP","Time","Rep")

plasmids<-cbind(variables,plasmids)

ggplot(plasmids,aes(x=VLP,y=value,color=Time))+
    geom_point()+geom_boxplot()+facet_grid(plasmid~.)
