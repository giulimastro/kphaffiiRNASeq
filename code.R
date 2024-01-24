#!/bin/python3
# grep "EMBL\tCDS" reference/genome.gff3 | cut -f 9 > genome.tab

# Change Folder Directory
directory<-"~/Github/kphaffii/"

setwd(directory)

homologs<-data.table::fread("homology_gs115/homology.txt",header=F)
keggtab <- data.table::fread("Pichia.blastkoala.txt",
                             header=T,
                             fill = T,
                             sep="\t",
                             na.strings="")



tab<-rtracklayer::readGFF("reference_cbs7435/CBS7435_21102016.gff3",
                          tags =c("ID",
                                  "Parent",
                                  "Dbxref",
                                  "Name",
                                  "gbkey",
                                  "gene",
                                  "locus_tag",
                                  "product",
                                  "protein_id",
                                  "inference") )
colnames(tab)<-c("ID",
                 "Parent",
                 "Dbxref",
                 "Name",
                 "gbkey",
                 "gene",
                 "locus_tag",
                 "product",
                 "protein_id",
                 "inference"
)

write.csv("annot.csv",
          tab,quote=F,
          row.names = F)

blast2go<-data.table::fread("blast2go_table.txt")

## RNAseq SK phaffii
# Said Muñoz Montero
# s.munoz-montero21@imperial.ac.uk

setwd("C:/Users/Giuli/OneDrive - Imperial College London/RNASeq/Analysis/kphaffii2/kphaffi_CHKV")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("DESeq2")
BiocManager::install("DESeq2")
library(DESeq2)

## Cite DESeq2
# Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq
# data with DESeq2 Genome Biology 15(12):550 (2014)

BiocManager::install("GO.db")
library(GO.db)
BiocManager::install("clusterProfiler")
library(clusterProfiler)

install.packages("tidyverse") #includes the ggplot2 & stringr
library(ggplot2)
library(stringr)
install.packages("pheatmap")
library(pheatmap)


# Pathway enrichment
BiocManager::install("enrichplot")
library(enrichplot)
BiocManager::install("pathview")
library(pathview)

BiocManager::install() #this installs all packages from BioConductor
.libPaths() #https://stackoverflow.com/questions/41839214/installation-path-not-writable-r-unable-to-update-packages
BiocManager::install(c("foreign", "KernSmooth", "Matrix", "mgcv", "nlme", "spatial", "survival"),
                     lib = "C:/Users/Giuli/AppData/Local/R/win-library/4.3")
#BiocManager::install("geneLenDataBase",lib = "C:/Users/Giuli/AppData/Local/R/win-library/4.3")
BiocManager::install("goseq",lib = "C:/Users/Giuli/AppData/Local/R/win-library/4.3")
library(goseq)
BiocManager::install("ballgown",lib = "C:/Users/Giuli/AppData/Local/R/win-library/4.3")
library(ballgown)

# Visualization
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
install.packages("UpSetR")
library(UpSetR)

# Formatting
install.packages("tidyr")
library(tidyr)
install.packages("RColorBrewer")
library(RColorBrewer)
install.packages("reshape2")
library("reshape2")
install.packages("dplyr")
library(dplyr)

install.packages("ggVennDiagram")
library(ggVennDiagram)
install.packages("circlize")
library(circlize)
BiocManager::install("AnnotationHub")
library(AnnotationHub)

install.packages("data.table", dependencies=TRUE)

# Loading counting files
#sampleFiles <- grep(".csv",list.files(directory),
#                    value=TRUE)

#sampleCondition=rep(c("CHKV18","GST","WT"),
#                    each=6)

#sampleTable <- data.frame(sampleName = sub(".csv","\\1",sampleFiles),
#                          fileName = sampleFiles,
#                          condition = sampleCondition)

#sampleTable$time<-factor(c(0,0,0,72,72,72))

# WT time

# sampleTableSubset <- subset(sampleTable, 
#                            condition %in% c("CHKV18"))

# ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTableSubset,
#                                       directory = directory,
#                                       design= ~ time )
#                                       design= ~ condition+time + condition:time)

# ddsHTSeq$time <- factor(ddsHTSeq$time)
# ddsHTSeq$time <- relevel(ddsHTSeq$time, 
#                          ref="0")

## DESeq normalization
#dds <- DESeq(ddsHTSeq, test="LRT", reduced= ~condition + time)

##### All Samples
ff <- list.files( path = "./counts", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
counts.files <- lapply( ff, read.table, skip = 4 )
# Column 4 is equivalent to htseq-count option -s reverse dUTP library prep
counts <- as.data.frame( sapply( counts.files, function(x) x[ , 4 ] ) )
ff <- gsub( "ReadsPerGene[.]out[.]tab", "", ff )
ff <- gsub( "[.]/counts/", "", ff )
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1

samples<-data.frame(Name=ff) %>% 
  separate(Name,c("VLP","Time","Rep"),"-",remove=T)

samples$condition <- paste(samples$VLP,samples$Time,sep="_")

rownames(samples)<-ff

homologs<-data.table::fread("homology_gs115/blast.txt",header=F)




# Empieza aqui :)

load("kphaffii.rdata")

blast2go<-blast2go[!duplicated(blast2go$Description),]

goTerms<-lapply(strsplit(blast2go$`GO IDs`,split="; "),substring,3,13)

IDs<-rep(blast2go$Description,unlist(lapply(goTerms, length)))

goseq.DF=data.frame(geneName=IDs,goTerms=unlist(goTerms))


load("counts.pichia.rdata")

load("counts_wCHKV18.rdata") #this one has the CHKV18 plasmid
#taking the plasmids from the counts file
plasmids<-data.frame(tail(counts,8)) #9 with the "counts_wCHKV18.rdata", otherwise its 8
plasmids$plasmid=rownames(plasmids)

plasmids<-melt(plasmids)
View(plasmids)

variables=str_split_fixed(plasmids$variable,"\\.",3)
colnames(variables)<-c("VLP","Time","Rep")

plasmids<-cbind(variables,plasmids)
View(plasmids)


ggplot(subset(plasmids,Time==72),aes(x=VLP,y=value,color=VLP))+
  geom_point()+facet_grid(plasmid~Time,scales="free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8,face="bold"))


ggplot(subset(plasmids, plasmid=="RABVG"),aes(x=VLP,y=value,color=VLP))+
  geom_point()+facet_grid(plasmid~Time,scales="free_y",) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8,face="bold"))

ggplot(plasmids,aes(x=VLP,y=value,color=VLP))+
  geom_point()+facet_grid(plasmid~Time,scales="free_y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8,face="bold"))

#Remove Rabies LZ and Rabies IV2
gotcha<-subset(subset(plasmids,VLP!="RabiesLZ"),VLP!="RabiesIV2")

#Rename RabiesFull to Rabies
#Rename CHKV18 to CHKV CAR
#Rename CHKV8 to CHKV ECSA
gotcha$VLP[gotcha$VLP == "CHKV18"] <- c("CHKV CAR")
gotcha$VLP[gotcha$VLP == "CHKV8"] <- c("CHKV ECSA")
gotcha$VLP[gotcha$VLP == "RabiesFull"] <- c("Rabies")

#Rename plasmid CHIKV08 to CHKV ECSA
gotcha$plasmid[gotcha$plasmid == "CHIKV08"] <- c("CHKV ECSA")
gotcha$plasmid[gotcha$plasmid == "CHIKV18"] <- c("CHKV CAR")


View(gotcha)

PlasmidsPlot<-ggplot(gotcha,aes(x=VLP,y=value,color=VLP))+
  geom_point()+facet_grid(plasmid~Time,scales="free_y") +
  ylab("Counts") +
  theme(strip.text.y.right = element_text(size=5)) +
  theme(axis.text.y = element_text(size=5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8,face="bold")) +
  theme(legend.position = "none")

PlasmidsPlot

svg("PlasmidsCountswithCHKVCAR.svg")
plot(PlasmidsPlot)
dev.off()


# Analysis

#If I want to remove the LZ and IV2 Rabies samples, do the following. Otherwise leave
#counts and samples as they are.

#Removing Rabies LZ and Rabies IV2
samples<-subset(subset(samples,VLP!="RabiesLZ"),VLP!="RabiesIV2")

#Rename RabiesFull to Rabies
#Rename CHKV18 to CHKV CAR
#Rename CHKV8 to CHKV ECSA
samples$VLP[samples$VLP == "CHKV18"] <- c("CHKV CAR")
samples$VLP[samples$VLP == "CHKV8"] <- c("CHKV ECSA")
samples$VLP[samples$VLP == "RabiesFull"] <- c("Rabies")


# Remove columns RabiesLZ and Rabies IV2


# Remove columns contains character
counts <-counts %>% select(-contains('RabiesLZ'))
counts <-counts %>% select(-contains('RabiesIV2'))
View(counts)

#This ends removing the Rabies LZ and IV2 from the dataset.

#Remove the plasmids from the counts for the Enhanced Volcano Plot
nrow(counts) ##5517 genes in total
tail(counts)
nrow(head(counts,-8)) #5509. Removing the plasmid sequences
tail(head(counts,-8)) #shows the last 6 genes in the df
counts1 <- head(counts, -8)
nrow(counts1) #verifying its 5509
tail(counts1) #verifying the last rows are not the plasmids.

#Do the DESeq2 analysis as follows:
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samples,
                              design= ~condition)
#Do the same but with counts1 fo the Enhanced Volcano plot
dds <- DESeqDataSetFromMatrix(countData = counts1,
                              colData = samples,
                              design= ~condition)
dds <- DESeq(dds)

vsd <- vst(dds, blind = F) #variance-stabilising transformation

#The returnData=True returns a dataframe to add labels to the points.
#If I wanted to add satisfying a certain condition see following link:
#https://stackoverflow.com/questions/15624656/label-points-in-geom-point

PCAplotDF<-plotPCA(vsd,intgroup=c("VLP","Time"),returnData=TRUE)
PCAplotDF$label<-NA
View(PCAplotDF)

PCAplotDF$label[PCAplotDF$name == "RabiesFull-72-A"] <- "Replicate A"
PCAplotDF$label[PCAplotDF$name == "RabiesFull-72-B"] <- "Replicate B"
PCAplotDF$label[PCAplotDF$name == "RabiesFull-72-C"] <- "Replicate C"

#if I don't want the labels on the graph, just remove the geom_text argument.

PCAplot<-plotPCA(vsd,intgroup=c("VLP","Time"))
PCAplot <- PCAplot +
  aes(color=VLP,shape=Time)+
  geom_text(data = PCAplotDF, aes(label = label),vjust=0, hjust=0)+
  theme_bw()+
  theme(legend.position = "bottom")+
  scale_color_manual(values=c(brewer.pal(3,"Dark2")[1:2],"blue",brewer.pal(8,"Dark2")[3:8],brewer.pal(3,"Paired")[1:2],"black"))

PCAplot


pdf("C:/Users/gpm15/OneDrive - Imperial College London/RNASeq/Analysis/PCA_RabiesFullOnly.pdf",width = 12)
plot(PCAplot)
dev.off()

svg("PCA_RabiesFullOnly_labelsRabies.svg",width = 12)
plot(PCAplot)
dev.off()


# Before norm
pdf("Boxplot.pdf")

boxplot(log2(counts(dds,normalize=F)+1),las=2, cex.axis=0.75, ylab="log2(counts)")

# After norm against transcript length
boxplot(log2(counts(dds,
                    normalize=T)+1),
        las=2,cex.axis=0.75, ylab="log2(norm counts)")
dev.off()

#Same, but in SVG
svg("C:/Users/gpm15/OneDrive - Imperial College London/RNASeq/Analysis/Boxplot.svg")
boxplot(log2(counts(dds,
                    normalize=F)+1),
        las=2, cex.axis=0.75, ylab="log2(counts)")
dev.off()

svg("C:/Users/gpm15/OneDrive - Imperial College London/RNASeq/Analysis/Boxplot_norm.svg")
# After norm against transcript length
boxplot(log2(counts(dds,
                    normalize=T)+1),
        las=2,cex.axis=0.75, ylab="log2(norm counts)")
dev.off()



resultsNames(dds)

##  Comparisons HPV 6, 11 & 16 vs GST & WT  ##

#72h
comparisonsCHKV<-data.frame("condition",A=c("HPV6_72","HPV6_72","HPV11_72","HPV11_72","HPV16_72","HPV16_72","GST_72"),B=c("GST_72","WT_72","GST_72","WT_72","GST_72","WT_72","WT_72"))
#0h
comparisonsCHKV<-data.frame("condition",A=c("HPV6_0","HPV6_0","HPV11_0","HPV11_0","HPV16_0","HPV16_0"),B=c("GST_0","WT_0","GST_0","WT_0","GST_0","WT_0"))
#HPVs against HPVs
comparisonsCHKV<-data.frame("condition",A=c("HPV6_0","HPV6_72","HPV6_0","HPV6_72","HPV16_0","HPV16_72"),B=c("HPV16_0","HPV16_72","HPV11_0","HPV11_72","HPV11_0","HPV11_72"))

View(comparisonsCHKV)

i=1
#Create a table with the first one and then append to that one
resCHKV <- results(dds,contrast = unlist(comparisonsCHKV[i,]))
head(resCHKV)
tail(resCHKV)
resCHKV@nrows
resCHKV@rownames


#selectFC1<- which(res$padj<0.05 & abs(res$log2FoldChange)>1)
selectFC2CHKV<- which(resCHKV$pvalue<0.05 & abs(resCHKV$log2FoldChange)>1)
#selectFC4<- which(res$padj<0.05 & abs(res$log2FoldChange)>4)
selectp005CHKV<- which(resCHKV$pvalue<0.05)
View(selectp005CHKV)
View(selectFC2CHKV)

resConditionCHKV<-resCHKV[selectFC2CHKV,]
#here I'm choosing the ones which have a p value smaller than 0.05 and an abs(log2FC)>1, meaning double the mRNA
head(resConditionCHKV)
nrow(resConditionCHKV)
View(resConditionCHKV)


resConditionCHKVp005<-resCHKV[selectp005CHKV,]
#choosing the ones with a p value smaller than 0.05
head(resConditionCHKVp005)
View(resConditionCHKVp005)
nrow(resConditionCHKVp005)


#Creating a list of the DEG of the second comparison so as to do the Venn Diagram
j=2
resCHKV2 <- results(dds,contrast = unlist(comparisonsCHKV[j,]))
selectFC2CHKV2<- which(resCHKV2$pvalue<0.05 & abs(resCHKV2$log2FoldChange)>1)
resConditionCHKV2<-resCHKV2[selectFC2CHKV2,]
View(resConditionCHKV2)
j=3
resCHKV3 <- results(dds,contrast = unlist(comparisonsCHKV[j,]))
selectFC2CHKV3<- which(resCHKV3$pvalue<0.05 & abs(resCHKV3$log2FoldChange)>1)
resConditionCHKV3<-resCHKV3[selectFC2CHKV3,]
View(resConditionCHKV3)
j=4
resCHKV4 <- results(dds,contrast = unlist(comparisonsCHKV[j,]))
selectFC2CHKV4<- which(resCHKV4$pvalue<0.05 & abs(resCHKV4$log2FoldChange)>1)
resConditionCHKV4<-resCHKV4[selectFC2CHKV4,]
View(resConditionCHKV4)
j=5
resCHKV5 <- results(dds,contrast = unlist(comparisonsCHKV[j,]))
selectFC2CHKV5<- which(resCHKV5$pvalue<0.05 & abs(resCHKV5$log2FoldChange)>1)
resConditionCHKV5<-resCHKV5[selectFC2CHKV5,]
View(resConditionCHKV5)
j=6
resCHKV6 <- results(dds,contrast = unlist(comparisonsCHKV[j,]))
selectFC2CHKV6<- which(resCHKV6$pvalue<0.05 & abs(resCHKV6$log2FoldChange)>1)
resConditionCHKV6<-resCHKV6[selectFC2CHKV6,]
View(resConditionCHKV6)

ListUpset<-list(HPV6vsGST_72h=resConditionCHKV@rownames,
        HPV6vsWT_72h=resConditionCHKV2@rownames,
        HPV11vsGST_72h=resConditionCHKV3@rownames,
        HPV11vsWT_72h=resConditionCHKV4@rownames,
        HPV16vsGST_72h=resConditionCHKV5@rownames,
        HPV16vsWT_72h=resConditionCHKV6@rownames)

ListUpsetHPV6<-list(HPV6vsGST_72h=resConditionCHKV@rownames,
                HPV6vsWT_72h=resConditionCHKV2@rownames)
ListUpsetHPV11<-list(HPV11vsGST_72h=resConditionCHKV3@rownames,
                    HPV11vsWT_72h=resConditionCHKV4@rownames)
ListUpsetHPV16<-list(HPV16vsGST_72h=resConditionCHKV5@rownames,
                    HPV16vsWT_72h=resConditionCHKV6@rownames)

#This venn diagram is for the resCHKV$pvalue<0.05 & abs(resCHKV$log2FoldChange)>1
#Basically: selectFC2CHKV index list

svg("VennDiagram_HPV6.svg")
ggVennDiagram(ListUpsetHPV6) +
  scale_fill_gradient(low="#F4FAFE",high = "#4981BF") +
  scale_color_brewer(palette = "Blues")
dev.off()

svg("VennDiagram_HPV11.svg")
ggVennDiagram(ListUpsetHPV11) +
  scale_fill_gradient(low="#F4FAFE",high = "#4981BF") +
  scale_color_brewer(palette = "Blues")
dev.off()

svg("VennDiagram_HPV16.svg")
ggVennDiagram(ListUpsetHPV16) +
  scale_fill_gradient(low="#F4FAFE",high = "#4981BF") +
  scale_color_brewer(palette = "Blues")
dev.off()

RColorBrewer::display.brewer.all()

#OR Upset
## Upset Gene level
CHKVUpset<-upset(fromList(ListUpset),order.by = "freq",nsets = 12)

svg("Upset_HPV.svg")
CHKVUpset
dev.off()

#to see the intersection, for example
intersect(ListUpset$CHKV8vsGST,ListUpset$CHKV8vsWT) #its different from 719 (from the upset graph) becaus eit also adds up all the other interesections as well.
reduce(intersect, ListUpset) #this gives the intersection of the 4 sub-lists

##VennDiagram & Upset of DEG finish here


##Heatmap of pvalue<0.05 
#https://github.com/mousepixels/sanbomics_scripts/blob/main/tutorial_complex_Heatmap.Rmd

df<-as.data.frame(resConditionCHKVp005)
head(df)
rownames(df)

df2<-as.data.frame(resConditionCHKV)
head(df2)
rownames(df2)

#The following table gives me the Gene name for each Feature, which is the rownames(df)
cbs7435.Annot<-data.table::fread("cbs7435_featurelist.csv", header=T)
View(cbs7435.Annot)

#Transform it into a list
keys<-cbs7435.Annot$Feature
values<-cbs7435.Annot$`Gene name`
keys
values
length(keys)
length(values)

l <- list()
for (i in 1:length(keys)){
  l[keys[i]] <-values[i]
}
l

#for non-mapped labels
no_values <- setdiff(rownames(df), keys)
for (i in 1:length(no_values)){
  l[no_values[i]] <- 'NA'
}
#adding this column of symbols
df$symbol <- unlist(l[rownames(df)], use.names = FALSE)
View(df)
nrow(df)

#take the top. The lower the baseMean, the noisier it will be. 50 is an arbitrary threshold.
#Play around with the log2FoldChange value, or it can be removed as well.

df.top <- df[(df$baseMean > 50) & (abs(df$log2FoldChange) > 0.5),]
df.top <- df[(abs(df$log2FoldChange) > 1),]
View(df.top)
nrow(df.top)

df.topALL <- df[ (df$baseMean > 50),]
View(df.topALL)
nrow(df.topALL)

#ordering
df.top <- df.top[order(df.top$log2FoldChange, decreasing = TRUE),]
View(df.top)
nrow(df.top)

df.topALL <- df.topALL[order(df.topALL$log2FoldChange, decreasing = TRUE),]
View(df.topALL)
nrow(df.topALL)

#removing the HEV, HPV18 and CHIKV08 from the data
rows_to_remove<-c("RABVG","CHIKV08","HPV18")
df.topALL<-df.topALL[!grepl(paste(rows_to_remove, collapse = "|"), rownames(df.topALL)), ]
nrow(df.topALL)

df.top<-df.top[!grepl(paste(rows_to_remove, collapse = "|"), rownames(df.top)), ]
nrow(df.top)


rlog_out <- rlog(dds, blind=FALSE) #get normalized count data from dds object
#faster alternative:
vsd <- vst(dds, blind = F)

#here its with df.top, but I can do it with df.topALL which has all the genes
mat<-assay(vsd)[rownames(df.top), colnames(vsd)] #sig genes x samples
colnames(mat) <- colnames(vsd)
View(mat)
colnames(mat)
nrow(mat)

#get the ones I want from colnames(mat)
#72h
keeps <- c("HPV11-72-A","HPV11-72-B","HPV11-72-C","WT-72-A","WT-72-B","WT-72-C","GST-72-A","GST-72-B","GST-72-C")
#0h
keeps<-c("HPV11-0-A","HPV11-0-B","HPV11-0-C","WT-0-A","WT-0-B","WT-0-C","GST-0-A","GST-0-B","GST-0-C")
mat2 <- as.data.frame(mat)[keeps]
View(mat2)
nrow(mat2)
as.matrix(mat2)
View(mat2)

base_mean <- rowMeans(mat2)
mat.scaled <- t(apply(mat2, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mat2)
View(mat.scaled)
nrow(mat.scaled)

#Getting the top 25 and the bottom 25 genes from the ordered df.top
#hence why it is so important to order it
num_keep <- 25
#1 to num_keep len-num_keep to len
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)) )

#Getting only the genes involved in WHAT i WANT
cellwall<-c("GDH3","DAS1","FLD","TDH3","FTR1","CZF1","PP7435_Chr2-0752")
cellwallRows<-df  %>%
  filter(symbol %in% cellwall)
rownames(cellwallRows) #gives me the row names of these ones, which are the genes
cellwallRows #this are the rows corresponding to the cellwall list of gene symbols
df[rownames(cellwallRows),] #another way to get the rows
cellwallRows$log2FoldChange #gives me the log2FC

#with the cell wall
l2_val <- as.matrix(cellwallRows$log2FoldChange)
colnames(l2_val)<-"logFC"
l2_val

#with the rows to keep
l2_val <- as.matrix(df.top[rows_keep,]$log2FoldChange) #getting log2 value for each gene we are keeping
colnames(l2_val)<-"logFC"

#with the cell wall 
mean <- as.matrix(cellwallRows$baseMean) #getting mean value for each gene we are keeping
colnames(mean)<-"AveExpr"
mean
#with the rows to keep
mean <- as.matrix(df.top[rows_keep,]$baseMean) #getting mean value for each gene we are keeping
colnames(mean)<-"AveExpr"



#maps values between b/w/r for min and max l2 values
col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red")) 

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                               height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F, 
              row_names_gp = gpar(fontsize = 3),
              column_names_gp = gpar(fontsize=4),
              column_labels = colnames(mat.scaled), name="Z-score",
              row_labels = df.top$symbol[rows_keep],
              cluster_columns = F)
h1 <- Heatmap(mat.scaled[rownames(cellwallRows),],
              row_labels = df.top[rownames(cellwallRows),]$symbol,
              row_names_gp = gpar(fontsize = 3),
              cluster_rows = F, 
              column_labels = colnames(mat.scaled), name="Z-score",
              cluster_columns = F) #changed the clustering of the columns to false
h1
h2 <- Heatmap(l2_val, row_labels = df.topALL[rownames(cellwallRows),]$symbol, 
              cluster_rows = F, name="logFC", top_annotation = ha, col = col_logFC,
              row_names_gp = gpar(fontsize = 3),
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(l2_val[i, j],2), x, y)
              })
h2 <- Heatmap(l2_val, row_labels = df.top$symbol[rows_keep], 
              cluster_rows = F, name="logFC", top_annotation = ha, col = col_logFC,
              row_names_gp = gpar(fontsize = 3),
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(l2_val[i, j],2), x, y, gp = gpar(fontsize=4))
              })
h2
h3 <- Heatmap(mean, row_labels = df.topALL[rownames(cellwallRows),]$symbol, 
              cluster_rows = F, name = "AveExpr", col=col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(mean[i, j],2), x, y)
              })
h3 <- Heatmap(mean, row_labels = df.top$symbol[rows_keep], 
              cluster_rows = F, name = "AveExpr", col=col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(mean[i, j],2), x, y)
              })
h3
h<-h1+h2#+h3
h

svg("Zscore_rabies_Time72.svg",width=12)
h1
dev.off()

#this is with the rows_keep
png("HEV_topDEG.png", res=250, width=1000, height=2000)
print(h)
dev.off()

#this is with the cell wall genes
png("HEV_cellwall.png", res=250, width=1000, height=2000)
print(h)
dev.off()

####


SummaryResults<-data.frame(ID=rownames(resConditionCHKV),
                           Up=sign(resConditionCHKV$log2FoldChange))
#Creating a data frame with ID and the sign of the log2fold change
colnames(SummaryResults)[2]<-paste(comparisonsCHKV[i,2],
                                   comparisonsCHKV[i,3],sep="-")
View(SummaryResults)
nrow(SummaryResults)

#to check the genes that are in p<0.05 but I want to see the log2fold change
FolChange_pCHKV<-data.frame(ID=rownames(resConditionCHKVp005),
                            log2FoldChange=resConditionCHKVp005$log2FoldChange)
colnames(FolChange_pCHKV)[2]<-paste("log2FoldChange",
                                    comparisonsCHKV[i,2],comparisonsCHKV[i,3],sep="-")
View(FolChange_pCHKV)
nrow(FolChange_pCHKV)

#removing HEV and HPV18 and CHIK08
#rows_to_remove<-c("HPV18","CHIKV08")
#FolChange_pCHKV<-FolChange_pCHKV[!grepl(paste(rows_to_remove, collapse = "|"), FolChange_pCHKV$ID), ]


#ordering
FolChange_pCHKV.ordered<-FolChange_pCHKV[order(FolChange_pCHKV$`log2FoldChange-HPV6_72-GST_72`,decreasing=T),]
View(FolChange_pCHKV.ordered)
nrow(FolChange_pCHKV.ordered)

resHomolog<-data.frame(GS115=homologs[match(x =rownames(resCHKV), homologs$V1),2],resCHKV)
#resExp <- assay(ntd)[selectFC2,]
View(resHomolog)
nrow(resHomolog)
resExp<-data.frame(resHomolog[selectFC2CHKV,])
View(resExp)
nrow(resExp)

## KEGG analysis
install.packages("remotes")
remotes::install_github("YuLab-SMU/GOSemSim")

#https://github.com/YuLab-SMU/clusterProfiler/issues/561
devtools::install_github("YuLab-SMU/clusterProfiler")
remotes::install_github("YuLab-SMU/createKEGGdb")
library(createKEGGdb)
species <-c("ppa")
createKEGGdb::create_kegg_db(species)
# You will get KEGG.db_1.0.tar.gz file in your working directory
install.packages("KEGG.db_1.0.tar.gz", repos=NULL,type="source")
library(KEGG.db)
KEGG.db::KEGG_dbfile()

kegg_gene_list<- resExp$log2FoldChange
names(kegg_gene_list)<-resExp$V2
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
head(kegg_gene_list)
names(kegg_gene_list)
length(kegg_gene_list)

keggEnrichment <- gseKEGG(geneList     = kegg_gene_list,
                          organism     = "ppa",
                          nPerm = 5500,
                          minGSSize    = 3,
                          maxGSSize    = 800,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "none",
                          use_internal_data = TRUE)

#removing the "- Komagataella phaffii" from each ID
keggEnrichment@result[,"Description"]<-sub("-.*$", "", keggEnrichment@result[,"Description"])
KEGGCHKV.df <- data.frame(keggEnrichment@result[,c("ID","Description","enrichmentScore")])
View(KEGGCHKV.df)
nrow(KEGGCHKV.df)
#result it's an object from keggEnrichment function
#the enrichment score is linked to the KEGG pathway. How much is that gene present in the pathway

colnames(KEGGCHKV.df)[3]<-paste(comparisonsCHKV[i,2],comparisonsCHKV[i,3],sep="-")
View(KEGGCHKV.df)


#I should do this within the loop, for each i, because the keggEnrichment changes.
KeggPlot<-dotplot(keggEnrichment, showCategory = 20, title = paste(comparisonsCHKV[i,2],comparisonsCHKV[i,3], sep="-"), split=".sign",font.size = 6) + facet_wrap(~.sign, ncol=1,scales = "free") 
#KeggPlot<-dotplot(keggEnrichment, showCategory = 20, title = "aaaa", split=".sign") + facet_wrap(~.sign, ncol=1,scales = "free") 

svg(paste0(paste("KEGG Enrichment",comparisonsCHKV[i,2],comparisonsCHKV[i,3], sep="-"),".svg"))
plot(KeggPlot)
#change y="pvalue" and the ylab="-Log10 P-value"
dev.off()

#KEGG pathview with the log2FC
pathview(gene.data  = kegg_gene_list,
                     pathway.id = "ppa00190",
                     species    = "ppa",
                     gene.idtype = "KEGG")
                     #limit      = list(gene=round(max(abs(kegg_gene_list))), cpd=1))


svg(filename=paste0("Volcano Plot ",(paste(comparisonsCHKV[i,2],comparisonsCHKV[i,3], sep="vs")),".svg"))
#change y="pvalue" and the ylab="-Log10 P-value"
EnhancedVolcano(resCHKV,
                lab=rownames(resCHKV),
                title=paste0("Volcano Plot ",(paste(comparisonsCHKV[i,2],comparisonsCHKV[i,3], sep="vs"))),
                x="log2FoldChange",
                y="padj",
                ylab = "-Log10 P-adj",
                #                selectLab = selectedGenes,
                pCutoff=0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 3.0,
                labCol= "black",
                labFace= "bold",
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = "none",
                drawConnectors = TRUE,
                caption="")
dev.off()

#This is for HPV6 vs WT, I don't change the i, but I do change the square brackets
resCHKV_EV <- results(dds,contrast = unlist(comparisonsCHKV[2,]))
svg(filename=paste0("Volcano Plot ",(paste(comparisonsCHKV[2,2],comparisonsCHKV[2,3], sep="vs")),".svg"))
EnhancedVolcano(resCHKV_EV,
                lab=rownames(resCHKV_EV),
                title=paste0("Volcano Plot ",(paste(comparisonsCHKV[2,2],comparisonsCHKV[2,3], sep="vs"))),
                x="log2FoldChange",
                y="padj",
                ylab = "-Log10 P-adj",
                #                selectLab = selectedGenes,
                pCutoff=0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 3.0,
                labCol= "black",
                labFace= "bold",
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = "none",
                drawConnectors = TRUE,
                caption="")
dev.off()

#if I want to check a specific gene, in this case PP7435_Chr2-0119
aaaa<-as.data.frame(resCHKV_EV)[as.data.frame(rownames(resCHKV_EV))==c("PP7435_Chr2-0119"),]
aaaa

#This is for HPV11 vs GST I don't change the i, but I do change the square brackets
resCHKV_EV <- results(dds,contrast = unlist(comparisonsCHKV[3,]))
svg(filename=paste0("Volcano Plot ",(paste(comparisonsCHKV[3,2],comparisonsCHKV[3,3], sep="vs")),".svg"))
EnhancedVolcano(resCHKV_EV,
                lab=rownames(resCHKV_EV),
                title=paste0("Volcano Plot ",(paste(comparisonsCHKV[3,2],comparisonsCHKV[3,3], sep="vs"))),
                x="log2FoldChange",
                y="padj",
                ylab = "-Log10 P-adj",
                #                selectLab = selectedGenes,
                pCutoff=0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 3.0,
                labCol= "black",
                labFace= "bold",
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = "none",
                drawConnectors = TRUE,
                caption="")
dev.off()

#This is for HPV11 vs WT, I don't change the i, but I do change the square brackets
resCHKV_EV <- results(dds,contrast = unlist(comparisonsCHKV[4,]))
svg(filename=paste0("Volcano Plot ",(paste(comparisonsCHKV[4,2],comparisonsCHKV[4,3], sep="vs")),".svg"))
EnhancedVolcano(resCHKV_EV,
                lab=rownames(resCHKV_EV),
                title=paste0("Volcano Plot ",(paste(comparisonsCHKV[4,2],comparisonsCHKV[4,3], sep="vs"))),
                x="log2FoldChange",
                y="padj",
                ylab = "-Log10 P-adj",
                #                selectLab = selectedGenes,
                pCutoff=0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 3.0,
                labCol= "black",
                labFace= "bold",
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = "none",
                drawConnectors = TRUE,
                caption="")
dev.off()

#This is for HPV16 vs GST, I don't change the i, but I do change the square brackets
resCHKV_EV <- results(dds,contrast = unlist(comparisonsCHKV[5,]))
svg(filename=paste0("Volcano Plot ",(paste(comparisonsCHKV[5,2],comparisonsCHKV[5,3], sep="vs")),".svg"))
EnhancedVolcano(resCHKV_EV,
                lab=rownames(resCHKV_EV),
                title=paste0("Volcano Plot ",(paste(comparisonsCHKV[5,2],comparisonsCHKV[5,3], sep="vs"))),
                x="log2FoldChange",
                y="padj",
                ylab = "-Log10 P-adj",
                #                selectLab = selectedGenes,
                pCutoff=0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 3.0,
                labCol= "black",
                labFace= "bold",
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = "none",
                drawConnectors = TRUE,
                caption="")
dev.off()

#This is for HPV16 vs WT, I don't change the i, but I do change the square brackets
resCHKV_EV <- results(dds,contrast = unlist(comparisonsCHKV[6,]))
svg(filename=paste0("Volcano Plot ",(paste(comparisonsCHKV[6,2],comparisonsCHKV[6,3], sep="vs")),".svg"))
EnhancedVolcano(resCHKV_EV,
                lab=rownames(resCHKV_EV),
                title=paste0("Volcano Plot ",(paste(comparisonsCHKV[6,2],comparisonsCHKV[6,3], sep="vs"))),
                x="log2FoldChange",
                y="padj",
                ylab = "-Log10 P-adj",
                #                selectLab = selectedGenes,
                pCutoff=0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 3.0,
                labCol= "black",
                labFace= "bold",
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = "none",
                drawConnectors = TRUE,
                caption="")
dev.off()


KEGGgenes.df<-list(data.frame(CBS7435=rownames(resExp),GS115=resExp$V2))
names(KEGGgenes.df)[1]<-paste("KEGG Genes",comparisonsCHKV[i,2],comparisonsCHKV[i,3],sep="-")

View(KEGGgenes.df)

## GO Analysis

# Upregulated genes
UpselectGO<- which(resCHKV$padj<0.05 & (resCHKV$log2FoldChange)>1)
#View(UpselectGO)
up.genes <- as.integer(blast2go$Description%in%rownames(resCHKV)[UpselectGO])
#View(up.genes)
names(up.genes)<-blast2go$Description
#up.genes<-up.genes[!(duplicated(names(up.genes)))]
up.goseq<-nullp(up.genes,bias.data = blast2go$Length)
#View(up.goseq)
up<-goseq(up.goseq,gene2cat = goseq.DF)
View(up)
#its strange that there are 2 terms in HPV16vsGST that dont have the term name. So I googled it
#and added it like this:
up[up$category=="GO:0016021",]$term <- "membrane"
up[up$category=="GO:0016021",]$ontology <- "CC"

up[upf$category=="GO:0031225",]$term <- "(obsolete) anchored component of membrane"
up[up$category=="GO:0031225",]$ontology <- "CC"



#setting the threshold of the pvalue to lower than 0.05
up<-subset(up,over_represented_pvalue<0.05)
upSample<-up[,c("category","term","ontology","over_represented_pvalue")]
colnames(upSample)[4]<-paste(comparisonsCHKV[i,2],comparisonsCHKV[i,3],sep="-")
#View(upSample)
upGO.df <- data.frame(upSample)
View(upGO.df)

#Visualising GO, top 15 hits
svg(paste0(paste("GO Up",comparisonsCHKV[i,2],comparisonsCHKV[i,3], sep="-"),".svg"),width=12)
up %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", title=paste(sub("_.*$", "", comparisonsCHKV[i,2])," vs ",sub("_.*$", "", comparisonsCHKV[i,3]),": GO Upregulated", sep=""), colour="p value", size="Count")
dev.off()

# Downregulated
DownselectGO<- which(resCHKV$padj<0.05 & (resCHKV$log2FoldChange)<1)
down.genes <- as.integer(blast2go$Description%in%rownames(resCHKV)[DownselectGO])
names(down.genes)<-blast2go$Description

#up.genes<-up.genes[!(duplicated(names(up.genes)))]

down.goseq<-nullp(down.genes,bias.data = blast2go$Length)
down<-goseq(down.goseq,gene2cat = goseq.DF)
down<-subset(down,over_represented_pvalue<0.05)
View(down)
downSample<-down[,c("category","term","ontology","over_represented_pvalue")]
colnames(downSample)[4]<-paste(comparisonsCHKV[i,2],comparisonsCHKV[i,3],sep="-")
downGO.df <- data.frame(downSample)
View(downGO.df)

svg(paste0(paste("GO Down",comparisonsCHKV[i,2],comparisonsCHKV[i,3], sep="-"),".svg"),width = 12)
down %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", title=paste(sub("_.*$", "", comparisonsCHKV[i,2])," vs ",sub("_.*$", "", comparisonsCHKV[i,3]),": GO Downregulated",sep=""), colour="p value", size="Count")
dev.off()

#here is where I change the i to the other rows because this for doesn't work.
for(i in 2:nrow(comparisonsCHKV)){
  resCHKV <- results(dds,contrast = unlist(comparisonsCHKV[i,]))
  
  #selectFC1<- which(res$padj<0.05 & abs(res$log2FoldChange)>1)
  selectFC2CHKV<- which(resCHKV$padj<0.05 & abs(resCHKV$log2FoldChange)>1)
  #selectFC4<- which(res$padj<0.05 & abs(res$log2FoldChange)>4)
  selectp005CHKV<- which(resCHKV$padj<0.05)
  
  resConditionCHKV<-resCHKV[selectFC2CHKV,]
  SummaryResults2<-data.frame(ID=rownames(resConditionCHKV), Up=sign(resConditionCHKV$log2FoldChange))
  colnames(SummaryResults2)[2]<-paste(comparisonsCHKV[i,2],comparisonsCHKV[i,3],sep="-")
  
  resConditionCHKVp005<-resCHKV[selectp005CHKV,]
  FolChange_pCHKV_2<-data.frame(ID=rownames(resConditionCHKVp005),log2FoldChange=resConditionCHKVp005$log2FoldChange)
  colnames(FolChange_pCHKV_2)[2]<-paste("log2FoldChange",comparisonsCHKV[i,2],comparisonsCHKV[i,3],sep="-")
  
  #Append
  SummaryResults<-merge(SummaryResults,SummaryResults2,all=T)
  FolChange_pCHKV<-merge(FolChange_pCHKV,FolChange_pCHKV_2,all=T)
  
  ## KEGG analysis append
  resHomolog<-data.frame(GS115=homologs[match(x =rownames(resCHKV), homologs$V1),2],resCHKV)
  # resExp <- assay(ntd)[selectFC2,]
  resExp <- data.frame(resHomolog[selectFC2CHKV,])
  
  ## KEGG analysis
  kegg_gene_list<- resExp$log2FoldChange
  names(kegg_gene_list)<-resExp$V2
  kegg_gene_list<-na.omit(kegg_gene_list)
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
  
  keggEnrichment <- gseKEGG(geneList     = kegg_gene_list,
                            organism     = "ppa",
                            nPerm = 5500,
                            minGSSize    = 3,
                            maxGSSize    = 800,
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "none",
                            use_internal_data = TRUE)
  
  #removing the "- Komagataella phaffii" from each ID
  keggEnrichment@result[,"Description"]<-sub("-.*$", "", keggEnrichment@result[,"Description"])
  
  KEGG.df2 <- data.frame(keggEnrichment@result[,c("ID","Description","enrichmentScore")])
  #View(KEGG.df2)
  colnames(KEGG.df2)[3]<-paste(comparisonsCHKV[i,2],comparisonsCHKV[i,3],sep="-")
  #View(KEGG.df2)
  KEGG.df2[is.na(KEGG.df2)]<-0
  KEGGCHKV.df<-merge(KEGGCHKV.df,KEGG.df2,all=T)
  #merging the original df with the one I'm creating here (df2)
  
  KEGGgenes.df2<-list(data.frame(CBS7435=rownames(resExp),GS115=resExp$V2))
  names(KEGGgenes.df2)[1]<-paste("KEGG Genes",comparisonsCHKV[i,2],comparisonsCHKV[i,3],sep="-")
  
  KEGGgenes.df<-c(KEGGgenes.df,KEGGgenes.df2)
  #merging the original df with the one I'm creating here (df2)
  
  #GO
  UpselectGO<- which(resCHKV$padj<0.05 & (resCHKV$log2FoldChange)>1)
  
  up.genes <- as.integer(blast2go$Description%in%rownames(resCHKV)[UpselectGO])
  
  names(up.genes)<-blast2go$Description
  
  #up.genes<-up.genes[!(duplicated(names(up.genes)))]
  
  up.goseq<-nullp(up.genes,bias.data = blast2go$Length)
  
  up<-goseq(up.goseq,gene2cat = goseq.DF)
  
  up<-subset(up,over_represented_pvalue<0.05)
  #View(up)
  upSample<-up[,c("category","term","ontology","over_represented_pvalue")]
  #View(upSample)
  colnames(upSample)[4]<-paste(comparisonsCHKV[i,2],comparisonsCHKV[i,3],sep="-")
  
  upGO.df<-merge(upGO.df,upSample,all=T)
  #View(upGO.df)
  # Down GO
  
  DownselectGO<- which(resCHKV$padj<0.05 & (resCHKV$log2FoldChange)<1)
  
  down.genes <- as.integer(blast2go$Description%in%rownames(resCHKV)[DownselectGO])
  
  names(down.genes)<-blast2go$Description
  
  #up.genes<-up.genes[!(duplicated(names(up.genes)))]
  
  down.goseq<-nullp(down.genes,bias.data = blast2go$Length)
  
  down<-goseq(down.goseq,gene2cat = goseq.DF)
  
  down<-subset(down,over_represented_pvalue<0.05)
  downSample<-down[,c("category","term","ontology","over_represented_pvalue")]
  colnames(downSample)[4]<-paste(comparisonsCHKV[i,2],comparisonsCHKV[i,3],sep="-")
  downGO.df<-merge(downGO.df,downSample,all=T)
  #View(downGO.df)
  
}
#everytime I finish running with each i, do the GO plot for up (949-960)
#and downregulated (977-988)
#and the KEGG plot (754-762 lines)

View(KEGGCHKV.df)
View(KEGGgenes.df)
View(downGO.df)
View(upGO.df)
View(SummaryResults)
View(FolChange_pCHKV)


upGO.df[is.na(upGO.df)]<-0
#colnames(upGO.df)[-1:-3]<-c("CHKV18vsGST","CHKV8vsGST", "CHKV8vsWT","CHKV18vsWT","CHKV8VSCHKV18")
View(upGO.df)
write.csv(upGO.df,file="UpGO_HPV6vs16vs11.csv")


downGO.df[is.na(downGO.df)]<-0
#colnames(downGO.df)[-1:-3]<-c("CHKV8vsGST","CHKV18vsGST", "CHKV8vsWT","CHKV18vsWT","CHKV8VSCHKV18")
View(downGO.df)
write.csv(downGO.df,file="DownGO_HPV6vs16vs11.csv")


#KEGG genes for each comparison
write.csv(KEGGgenes.df$`KEGG Genes-HPV6_72-GST_72`,file="KEGGGenes_HPV6vGST.csv")
write.csv(KEGGgenes.df$`KEGG Genes-HPV6_72-WT_72`,file="KEGGGenes_HPV6vWT.csv")
write.csv(KEGGgenes.df$`KEGG Genes-HPV11_72-GST_72`,file="KEGGGenes_HPV11vGST.csv")
write.csv(KEGGgenes.df$`KEGG Genes-HPV11_72-WT_72`,file="KEGGGenes_HPV11vWT.csv")
write.csv(KEGGgenes.df$`KEGG Genes-HPV16_72-GST_72`,file="KEGGGenes_HPV16vGST.csv")
write.csv(KEGGgenes.df$`KEGG Genes-HPV16_72-WT_72`,file="KEGGGenes_HPV16vWT.csv")


FolChange_pCHKV[is.na(FolChange_pCHKV)]<-0
FolChange_pCHKV.ordered<-FolChange_pCHKV[order(FolChange_pCHKV$`log2FoldChange-HPV6_72-HPV16_72`,decreasing=T),]
View(FolChange_pCHKV.ordered)
write.csv(FolChange_pCHKV.ordered,file="log2FoldChange_HPV6vs16vs11.csv")


#For the Metabolic Model
mmodel<-data.table::fread("Modelgenes.csv", header=T)
View(mmodel[,2])

View(FolChange_pCHKV)
View(homologs)
FolChange_pCHKVHom<-FolChange_pCHKV
View(FolChange_pCHKVHom)

#This is an example
#FolChange_pCHKV$ID[2628] the same as: FolChange_pCHKV[2628,1]
FolChange_pCHKVHom$ID[2628]<-homologs$V2[match(x=FolChange_pCHKVHom$ID[2628],homologs$V1)]
FolChange_pCHKVHom[2628,]
View(FolChange_pCHKVHom)

p=1
FolChange_pCHKVHom$ID[p]<-homologs$V2[match(x=FolChange_pCHKVHom$ID[p],homologs$V1)]

for (p in 2:nrow(FolChange_pCHKVHom)){
  FolChange_pCHKVHom$ID[p]<-homologs$V2[match(x=FolChange_pCHKVHom$ID[p],homologs$V1)]
}
FolChange_pCHKVHom<-na.omit(FolChange_pCHKVHom) #removing the NAs (the homologs that don't exist)
View(FolChange_pCHKVHom)

m=1
IDmmodel<-FolChange_pCHKVHom[match(x=mmodel$`0`[m],FolChange_pCHKVHom$ID),]
View(IDmmodel)

for (m in 2:nrow(mmodel)){
  IDmmodel2<-FolChange_pCHKVHom[match(x=mmodel$`0`[m],FolChange_pCHKVHom$ID),]
  IDmmodel<-merge(IDmmodel,IDmmodel2,all=T)
}
write.csv(IDmmodel,file="IDsForModel_HEV.csv")



#Read the Annotated pichia genome
cbs7435.Annot<-data.table::fread("cbs7435_featurelist.csv", header=T)
View(cbs7435.Annot)

### Forma para ver si hay ALGUN gen en particular que tenga p<0.05, te devuelve la posicion y los valores
match(x="PP7435_Chr2-0119",FolChange_pCHKV.ordered$ID)
FolChange_pCHKV.ordered[match(x="PP7435_Chr2-0119",FolChange_pCHKV.ordered$ID),1:ncol(FolChange_pCHKV.ordered)]

match(x="PP7435_Chr2-0119",cbs7435.Annot$Feature)
cbs7435.Annot[match(x="PP7435_Chr2-0119",cbs7435.Annot$Feature),1:ncol(cbs7435.Annot)]
cbs7435.Annot$`Gene description`[match(x="PP7435_Chr2-0119",cbs7435.Annot$Feature)]
### END ###

FolChange_pCHKV.ordered$ID[1]
cbs7435.Annot$`Gene description`[match(x=FolChange_pCHKV.ordered$ID[1],cbs7435.Annot$Feature)]

#Creating a dataframe with the log2Fold change and including the Annotation.
j=1
AnnotLog2FC<-data.frame(ID=FolChange_pCHKV.ordered$ID[j],
           Gene.Description=cbs7435.Annot$`Gene description`[match(x=FolChange_pCHKV.ordered$ID[j],cbs7435.Annot$Feature)],
           Gene.name=cbs7435.Annot$`Gene name`[match(x=FolChange_pCHKV.ordered$ID[j],cbs7435.Annot$Feature)],
           log2FC_HPV6vsGST=FolChange_pCHKV.ordered$`log2FoldChange-HPV6_72-GST_72`[j],
           log2FC_HPV6vsWT=FolChange_pCHKV.ordered$`log2FoldChange-HPV6_72-WT_72`[j],
           log2FC_HPV11vsGST=FolChange_pCHKV.ordered$`log2FoldChange-HPV11_72-GST_72`[j],
           log2FC_HPV11vsWT=FolChange_pCHKV.ordered$`log2FoldChange-HPV11_72-WT_72`[j],
           log2FC_HPV16vsGST=FolChange_pCHKV.ordered$`log2FoldChange-HPV16_72-GST_72`[j],
           log2FC_HPV16vsWT=FolChange_pCHKV.ordered$`log2FoldChange-HPV16_72-WT_72`[j]
           )
View(AnnotLog2FC)

#for HPV6 vs HPV16 vs HPV11
j=1
AnnotLog2FC<-data.frame(ID=FolChange_pCHKV.ordered$ID[j],
                        Gene.Description=cbs7435.Annot$`Gene description`[match(x=FolChange_pCHKV.ordered$ID[j],cbs7435.Annot$Feature)],
                        Gene.name=cbs7435.Annot$`Gene name`[match(x=FolChange_pCHKV.ordered$ID[j],cbs7435.Annot$Feature)],
                        log2FC_HPV6vsHPV16_Time0=FolChange_pCHKV.ordered$`log2FoldChange-HPV6_0-HPV16_0`[j],
                        log2FC_HPV6vsHPV16=FolChange_pCHKV.ordered$`log2FoldChange-HPV6_72-HPV16_72`[j],
                        log2FC_HPV6vsHPV11_Time0=FolChange_pCHKV.ordered$`log2FoldChange-HPV6_0-HPV11_0`[j],
                        log2FC_HPV6vsHPV11=FolChange_pCHKV.ordered$`log2FoldChange-HPV6_72-HPV11_72`[j],
                        log2FC_HPV16vsHPV11_Time0=FolChange_pCHKV.ordered$`log2FoldChange-HPV16_0-HPV11_0`[j],
                        log2FC_HPV16vsHPV11=FolChange_pCHKV.ordered$`log2FoldChange-HPV16_72-HPV11_72`[j]
)
View(AnnotLog2FC)


for (j in 2:nrow(FolChange_pCHKV.ordered)){
  AnnotLog2FC.2<-data.frame(ID=FolChange_pCHKV.ordered$ID[j],
                             Gene.Description=cbs7435.Annot$`Gene description`[match(x=FolChange_pCHKV.ordered$ID[j],cbs7435.Annot$Feature)],
                             Gene.name=cbs7435.Annot$`Gene name`[match(x=FolChange_pCHKV.ordered$ID[j],cbs7435.Annot$Feature)],
                             log2FC_HPV6vsGST=FolChange_pCHKV.ordered$`log2FoldChange-HPV6_72-GST_72`[j],
                             log2FC_HPV6vsWT=FolChange_pCHKV.ordered$`log2FoldChange-HPV6_72-WT_72`[j],
                             log2FC_HPV11vsGST=FolChange_pCHKV.ordered$`log2FoldChange-HPV11_72-GST_72`[j],
                             log2FC_HPV11vsWT=FolChange_pCHKV.ordered$`log2FoldChange-HPV11_72-WT_72`[j],
                             log2FC_HPV16vsGST=FolChange_pCHKV.ordered$`log2FoldChange-HPV16_72-GST_72`[j],
                             log2FC_HPV16vsWT=FolChange_pCHKV.ordered$`log2FoldChange-HPV16_72-WT_72`[j]
                            )
  AnnotLog2FC<-merge(AnnotLog2FC,AnnotLog2FC.2,all=T)
}

#for HPV6 vs HPV16 vs HPV11
for (j in 2:nrow(FolChange_pCHKV.ordered)){
  AnnotLog2FC.2<-data.frame(ID=FolChange_pCHKV.ordered$ID[j],
                            Gene.Description=cbs7435.Annot$`Gene description`[match(x=FolChange_pCHKV.ordered$ID[j],cbs7435.Annot$Feature)],
                            Gene.name=cbs7435.Annot$`Gene name`[match(x=FolChange_pCHKV.ordered$ID[j],cbs7435.Annot$Feature)],
                            log2FC_HPV6vsHPV16_Time0=FolChange_pCHKV.ordered$`log2FoldChange-HPV6_0-HPV16_0`[j],
                            log2FC_HPV6vsHPV16=FolChange_pCHKV.ordered$`log2FoldChange-HPV6_72-HPV16_72`[j],
                            log2FC_HPV6vsHPV11_Time0=FolChange_pCHKV.ordered$`log2FoldChange-HPV6_0-HPV11_0`[j],
                            log2FC_HPV6vsHPV11=FolChange_pCHKV.ordered$`log2FoldChange-HPV6_72-HPV11_72`[j],
                            log2FC_HPV16vsHPV11_Time0=FolChange_pCHKV.ordered$`log2FoldChange-HPV16_0-HPV11_0`[j],
                            log2FC_HPV16vsHPV11=FolChange_pCHKV.ordered$`log2FoldChange-HPV16_72-HPV11_72`[j]
  )
  AnnotLog2FC<-merge(AnnotLog2FC,AnnotLog2FC.2,all=T)
}

#AnnotLog2FC<-AnnotLog2FC[,-ncol(AnnotLog2FC)] #removing the GEne.Description1 column at the end
View(AnnotLog2FC)
#changing the name of the column to symbol
colnames(AnnotLog2FC)[colnames(AnnotLog2FC) == "Gene.name"] <- "symbol"

write.csv(AnnotLog2FC,file="log2FoldChange_HPV6vs16vs11_Annot.csv")

##Creating a heatmap of the log2FC##

#loading the csv file so as not to do EVERYTHING
AnnotLog2FC<-read.csv(file = "log2FoldChange_HPV_Annot.csv")
View(AnnotLog2FC)

#Getting only the genes involved in what I'm interested in
cellwall<-c("GDH3","DAS1","FLD","TDH3","FTR1","CZF1","PP7435_Chr2-0752")
cellwall<-c("TRR1","YAP1","PMP20","SOD1","PRX1-2","YPD1","HYR1","AHP1","SKN7","PRX1-1","PRX5") #oxidative stress
cellwall<-c("RPT6","RPN1","RPN2","RPN8","RPN11","RPN9","RPN12","RPN13","RPN14","UBX2","NPL4","UBX4","CDC48","UBC1","PP7435_Chr1-0354","HRD3","DER1","DOA1") #proteasome
cellwall<-c("LHS1","JEM1","KAR2","PDI1","IRE1","CNE1")#UPR markers
cellwall<-c("PEX1","PEX11","PEX11C","PEX10","PEX12","PEX13","PEX14","PEX17","PEX19","PEX2","PEX20","PEX22","PEX25","PEX28","PEX29","PEX3","PEX30","PEX31","PEX4","PEX5","PEX6","PEX7","PEX8","PEX26")#peroxisome "AOX2","AOX1","CTA1","DAK2","DAS1","DAS2","FDH1","FGH1","FLD"
cellwall<-c("YIH1","GCD11","CDC123","PRT1","HCR1","NIP1","TIF35","TIF34","RLI1","TIF2","TIF4632") #TRANSLATION
cellwall<-c("UBX2","CDC48","UFD1","SSM4","HRD3","PP7435_Chr1-0354","CUZ1","DOA1","VMS1","UBX4","HRD1-2","HRD1-1","UBX5","RPT4","HLJ1","UBC1","EPS1","SBH1") #ERAD


cellwallRows<-AnnotLog2FC  %>%
  filter(symbol %in% cellwall)
View(cellwallRows)
names(cellwallRows)[names(cellwallRows) == "log2FC_HPV11vsGST"] <- "HPV-11 vs GST" #changing the name of the rows
names(cellwallRows)[names(cellwallRows) == "log2FC_HPV11vsWT"] <- "HPV-11 vs WT"


rownames(cellwallRows)<- cellwallRows$X #putting the X as the rowname
cellwallRows<-cellwallRows[,-1] #i need this when I load the AnnotLog2FC csv file

names(cellwallRows)
cellwallmatrix<-as.matrix(cellwallRows[,6:7])
cellwallmatrix<-as.matrix(cellwallRows[,c("HPV-11 vs WT","HPV-11 vs GST")])
View(cellwallmatrix)

ComplexHeatmap::pheatmap(cellwallmatrix,
                         labels_row = cellwallRows$symbol,
                         fontsize = 5,
                         cluster_rows = F,
                         cluster_cols = F,
                         heatmap_legend_param = list(title = "log2FC", at = c(round(min(cellwallmatrix),2),round(max(cellwallmatrix),2)))
                         )


#Annotated Genes with p<0.05 and abs(log2FC)>1 with their annotation
#These are the ones used for KEGG (resCHKV)
AnnotGenesKEGG<-subset(AnnotLog2FC,abs(AnnotLog2FC[,4:ncol(AnnotLog2FC)])>1)
View(AnnotGenesKEGG)
nrow(AnnotGenesKEGG)
nrow(AnnotLog2FC)

# KEGG pathways in excel format AND heatmap
KEGGCHKV.df[is.na(KEGGCHKV.df)]<-0
View(KEGGCHKV.df)
#KEGGgenes<-data.frame(rownames(resCHKV))
colnames(KEGGCHKV.df)[-1:-2]<-c("HPV6vsGST","HPV6vsWT","HPV11vsGST","HPV11vsWT","HPV16vsGST","HPV16vsWT")
View(KEGGCHKV.df)
write.csv(KEGGCHKV.df,file="KEGGHPV.csv")

#KEGGCHKV.df<-KEGGCHKV.df[,c("ID","Description","RabiesvsGST","RabiesvsWT","RabiesvsGST_Time0","RabiesvsWT_Time0")]
pdf("KEGGHPV_heatmap.pdf",height=20)
svg("KEGGHPV_heatmap.svg",width=12)
pheatmap(KEGGCHKV.df[,-1:-2],
         labels_row = paste(KEGGCHKV.df$ID,KEGGCHKV.df$Description, sep=" - "),
         #color=c("#C8102E","white","#006341"),
         fontsize=5,
         cluster_cols = F)
dev.off()

#same heatmap but only 72h, removing the ppa01240 because this one is only in the Time 0
KEGGCHKV.df.HM<-KEGGCHKV.df[KEGGCHKV.df$ID!="ppa01240",]
View(KEGGCHKV.df.HM)
pheatmap(KEGGCHKV.df.HM[,3:4],
         labels_row = paste(KEGGCHKV.df.HM$ID,KEGGCHKV.df.HM$Description, sep=" - "),
         #color=c("#C8102E","white","#006341"),
         fontsize=5,
         cluster_cols = FALSE)
dev.off()

#Another visualisation of the KEGG (ERROR)
pdf("KEGGPlots.pdf")
for (a in 1:4){
  KeggPlot<-dotplot(keggEnrichment, showCategory = 20, title = comparisonsCHKV[a,3] , split=".sign") + facet_wrap(~.sign, ncol=1,scales = "free") 
  plot(KeggPlot)
}
dev.off()

KeggPlot<-dotplot(keggEnrichment, showCategory = 20, title = paste0(comparisonsCHKV[i,2],comparisonsCHKV[i,3], sep="-") , split=".sign") + facet_wrap(~.sign, ncol=1,scales = "free") 
pdf(paste(comparisonsCHKV[i,2],comparisonsCHKV[i,3],".pdf", sep= "-"))
plot(KeggPlot)
dev.off()


SummaryResults[is.na(SummaryResults)]<-0
View(SummaryResults)

rownames(SummaryResults)<-SummaryResults$ID

SummaryResults<-SummaryResults[,-1] #taking out the ID column
View(SummaryResults)
colnames(SummaryResults)<-c("GSTvsCHKV8","GSTvsCHKV18","GSTvsWT","CHKV8vsWT")
colnames(SummaryResults)<-c("HEVvsWT", "HEVvsGST")

Result<-SummaryResults
Result.ordered<-Result[order(Result$HEVvsWT,decreasing = T),]
View(Result)
View(Result.ordered)

rownames(Result.ordered)
cbs7435.Annot$`Gene description`[match(x=rownames(Result.ordered)[1],cbs7435.Annot$Feature)]

#Creating a dataframe with the Result ordered with the genes from with p<0.05 and abs(log2FC)>1 & including the Annotation ONLY for WT
k=1
AnnotResult.ordered<-data.frame(ID=rownames(Result.ordered)[k],
                        Gene.Description=cbs7435.Annot$`Gene description`[match(x=rownames(Result.ordered)[k],cbs7435.Annot$Feature)])
View(AnnotResult.ordered)

for (k in 2:nrow(Result.ordered)){
  AnnotResult.ordered2<-data.frame(ID=rownames(Result.ordered)[k],
                            Gene.Description=cbs7435.Annot$`Gene description`[match(x=rownames(Result.ordered)[k],cbs7435.Annot$Feature)],
                            Gene.name=cbs7435.Annot$`Gene name`[match(x=rownames(Result.ordered)[k],cbs7435.Annot$Feature)])
  AnnotResult.ordered<-merge(AnnotResult.ordered,AnnotResult.ordered2,all=T)
}
View(AnnotResult.ordered)

rownames(AnnotResult.ordered)<-AnnotResult.ordered$ID
Result.ordered.Annot<-merge(Result.ordered,AnnotResult.ordered,by=0, all=T)
View(Result.ordered.Annot)
rownames(Result.ordered.Annot)<-Result.ordered.Annot$ID
Result.ordered.Annot<-Result.ordered.Annot[,-1]
Result.ordered.Annot<-Result.ordered.Annot[,c(4,6,5,1,2,3)] #reordering columns
View(Result.ordered.Annot)
write.csv(Result.ordered.Annot,file="UpAndDownGenesCHKV.csv")
##IT ISNT ORDERED. Check


#Forma para ver si hay ALGUN gen en particular que tenga p<0.05 and abs(log2FoldChange)>1, te devuelve la posicion y los valores
match(x="PP7435_Chr2-0001",SummaryResults$ID)
View(SummaryResults[match(x="PP7435_Chr2-0001",SummaryResults$ID),1:ncol(SummaryResults)])

#Contrast<-Result-Result$WT

pdf("heatmap_GSTvsCHKV8_18_WT.pdf",width = 12,height = 9)
#pheatmap(Result-Result$WT,cluster_rows = T,show_rownames = F,legend_breaks = c(-2,0,2),color=c("#006341","white","#C8102E"))
pheatmap(Result.ordered,cluster_rows = F,
         show_rownames = F,legend_breaks = c(-1,0,1),
         color=c("#006341","white","#C8102E"))

dev.off()

## Nivel Gen from KEGG analysis
genes=rownames(Result)
#taking the gene names from the Result data frame
View(genes)
View(Result)

#ForHEV
Up<-list(HEVvsWT=genes[Result$HEVvsWT>0],
         HEVvsGST=genes[Result$HEVvsGST>0])
#ForRabies
Up<-list(WT=genes[Result$WT>0],
         GST=genes[Result$GST>0],
         RabiesIV2=genes[Result$RabiesIV2>0],
         RabiesLZ=genes[Result$RabiesLZ>0])
View(Up)

#ForHEV
Down<-list(HEVvsWT=genes[Result$HEVvsWT<0],
         HEVvsGST=genes[Result$HEVvsGST<0])
#ForRabies
Down<-list(GST=genes[Result$GST<0],
           RabiesIV2=genes[Result$RabiesIV2<0],
           RabiesLZ=genes[Result$RabiesLZ<0],
           WT=genes[Result$WT<0])
View(Down)

pdf("upset.pdf")
upset(fromList(Up),order.by = "freq",nsets = 12)
upset(fromList(Down),order.by = "freq",nsets = 12)
dev.off()

## Nivel KEGG Pathways
pathways=KEGGRabies.df$ID
#taking the pathways ID from the KEGG analysis

UpPathways<-list(GST=pathways[KEGGRabies.df$GST>0],
                RabiesIV2=pathways[KEGGRabies.df$RabiesIV2>0],
                 RabiesLZ=pathways[KEGGRabies.df$RabiesLZ>0],
                 WT=pathways[KEGGRabies.df$WT>0])
View(UpPathways)

UpPath_WT<-data.frame(UpPathways$WT)
UpPath_GST<-data.frame(UpPathways$GST)
UpPath_LZ<-data.frame(UpPathways$RabiesLZ)
UpPath_IV2<-data.frame(UpPathways$RabiesIV2)


write.table(UpPath_WT,file="UpPathRabvsWT_GST_LZ_IV2.csv")
write.table(UpPath_GST,"UpPathRabvsWT_GST_LZ_IV2.csv",append=TRUE)
write.table(UpPath_LZ,"UpPathRabvsWT_GST_LZ_IV2.csv",append=TRUE)
write.table(UpPath_IV2,"UpPathRabvsWT_GST_LZ_IV2.csv",append=TRUE)



DownPathways<-list(GST=pathways[KEGGRabies.df$GST<0],
                   RabiesIV2=pathways[KEGGRabies.df$RabiesIV2<0],
                   RabiesLZ=pathways[KEGGRabies.df$RabiesLZ<0],
                   WT=pathways[KEGGRabies.df$WT<0])
View(DownPathways)

DownPath_WT<-data.frame(DownPathways$WT)
View(Down_WT)
DownPath_GST<-data.frame(DownPathways$GST)
DownPath_LZ<-data.frame(DownPathways$RabiesLZ)
DownPath_IV2<-data.frame(DownPathways$RabiesIV2)


write.table(DownPath_WT,file="DownPathRabvsWT_GST_LZ_IV2.csv")
write.table(DownPath_GST,"DownPathRabvsWT_GST_LZ_IV2.csv",append=TRUE)
write.table(DownPath_LZ,"DownPathRabvsWT_GST_LZ_IV2.csv",append=TRUE)
write.table(DownPath_IV2,"DownPathRabvsWT_GST_LZ_IV2.csv",append=TRUE)


pdf("upset.pdf")
upset(fromList(UpPathways),order.by = "freq",nsets = 12)
upset(fromList(DownPathways),order.by = "freq",nsets = 12)
dev.off()

## Nivel GO Enrichment
GOTermDown=downGO.df$term
GOTermUp=upGO.df$term

#taking the pathways term from the GO Enrichment

UpEnrichment<-list(HEVvsWT=GOTermUp[upGO.df$HEVvsWT>0],
                   HEVvsGST=GOTermUp[upGO.df$HEVvsGST>0])
DownEnrichment<-list(HEVvsWT=GOTermDown[downGO.df$HEV_72.WT_72>0],
                   HEVvsGST=GOTermDown[downGO.df$`HEV_72-GST_72`>0])

View(UpEnrichment)
View(DownEnrichment)

pdf("upsetGO_HEV.pdf")
upset(fromList(UpEnrichment),order.by = "freq",nsets = 12)
upset(fromList(DownEnrichment),order.by = "freq",nsets = 12)
#in the case of HEV, given there are no downregulated enriched pathways, this upset can't be done
dev.off()


#This is the volcano plot for the i=2 (should do this for all i)
svg(paste0(paste("Volcano Plot ",sub("_.*$", "", comparisonsCHKV[i,2])," vs ",sub("_.*$", "", comparisonsCHKV[i,3]), sep=""),".svg"))
EnhancedVolcano(resCHKV,
                lab=rownames(resCHKV),
                title=paste("Volcano Plot ",sub("_.*$", "", comparisonsCHKV[i,2])," vs ",sub("_.*$", "", comparisonsCHKV[i,3]), sep=""),
                #title=paste("Volcano Plot",comparisonsCHKV[i,2],comparisonsCHKV[i,3], sep="-"),
                x="log2FoldChange",
                y="padj",
                ylab = "-Log10 P-adj",
                #                selectLab = selectedGenes,
                pCutoff=0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 3.0,
                labCol= "black",
                labFace= "bold",
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = "none",
                drawConnectors = TRUE,
                caption="")
dev.off()

#Cant use enrichGO from clustalProfiler because it doesnt have pichia's database
go_enrich<- enrichGO(gene = kegg_gene_list,
                     universe = names(kegg_gene_list),
                     OrgDb = 
                     #keyType = 'ENSEMBL',
                     readable = T,
                     ont = "ALL",
                     pvalueCutoff = 0.05, 
                     qvalueCutoff = 0.10)

go_enrich

ah <- AnnotationHub()
query(ah, "OrgDb")
query(ah, c("OrgDb","pichia"))
query(ah, "pichia")
#Doesn't have Orgdb, it has inparanoid8db which is a database of orthologs and paralogs
