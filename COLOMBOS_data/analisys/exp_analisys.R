install.packages("BiocManager")
install.packages("forcats")
install.packages("stringr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readr")
install.packages("tidyr")
install.packages("survminer")
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("org.Hs.eg.db")



#qunatile
BiocManager::install("preprocessCore")
#Para usar ComBat. Los instale from source
BiocManager::install("genefilter")
BiocManager::install("sva")

# Para rma
BiocManager::install("metafor")

# Para normalize
BiocManager::install("EnrichmentBrowser")



#Recursos
#https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html#Filtering_lowly_expressed_genes
#https://www.stat.purdue.edu/bigtap/online/docs/Introduction_to_Microarray_Analysis_GSE15947.html
#https://sbc.shef.ac.uk/geo_tutorial/tutorial.nb.html
#https://rpubs.com/sharbie/462299

#GEOquery
####### 
library(GEOquery)
library(limma)
## change my_id to be the dataset that you want.
my_id <- "GSE1642"
gse <- getGEO(my_id)

length(gse)
gse_GLP1343 <- gse[[1]]
gse_GLP1396<- gse[[2]]

boxplot(exprs(gse_GLP1343),outline=FALSE)


#######
# Manual
# exp_PROM_Original.txt
expTable_PROM_original <- read.table("exp_PROM_Original.txt", header=TRUE, row.names=1, sep=",")

# exp colombos GSE1642
colombos_GSE1642 <- read.table("/home/agustin/FBA_Tesis/PROM_trabajo/COLOMBOS_data/analisys/colombos_GSE1642.txt",  row.names=1, dec = "." ,skip=7,  header=FALSE, fill=TRUE)
colombos_GSE16811 <- read.table("/home/agustin/FBA_Tesis/PROM_trabajo/COLOMBOS_data/analisys/colombos_GSE16811.txt",  row.names=1, dec = "." ,skip=7,  header=FALSE, fill=TRUE)
colombos_9776 <- read.table("/home/agustin/FBA_Tesis/PROM_trabajo/COLOMBOS_data/colombos_9776.txt",  row.names=1, dec = "." ,skip=7,  header=FALSE, fill=TRUE)
colombos_16811_40846_9776_35362_14840 <- read.table("/home/agustin/FBA_Tesis/PROM_trabajo/COLOMBOS_data/colombos_16811_40846_9776_35362_14840.txt",  row.names=1, dec = "." ,skip=7,  header=FALSE, fill=TRUE)

select <- as.matrix (colombos_GSE1642)
boxplot(select,outline=FALSE)

# limma 
y2 <- normalizeBetweenArrays(log2(select), method="quantile")
boxplot(y2,outline=FALSE)

#quantile
library("preprocessCore")
select_quantile=normalize.quantiles(select)
boxplot(select_quantile,outline=FALSE)

# # rma
# library("metafor")
# data.rma = rma(colombos_GSM27855)
# data.matrix = exprs(data.rma)
# 
# # normalize
# library("EnrichmentBrowser")
# maSE <- normalize(colombos_GSM27855) 

#ComBat

a=replicate(138, 1)
b=replicate(299, 2)
# c=replicate(17, 3)
# d=replicate(16, 4)
# e=replicate(6,5)
batch<- c(a,b)
length(batch)
dim(select_quantile)

library("sva")
select_quantile[is.na(select_quantile)] = 0
data_combat=ComBat(dat=select_quantile, batch=batch, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
boxplot(data_combat,outline=FALSE)

#####
# MDS
mds = plotMDS(expTable_PROM_original,xlim=c(-2,2) ,xlab = "Dim1", ylim=c(-2,2), ylab = "Dim2", cex.lab= 1.5,pch = 15,cex=2, col=c( rep("black",6),  rep("red",12), rep("yellow",12), rep("orange",12))) #,rep("green",3),rep("deeppink",3),rep("gray",3), rep("pink",3)))
plotMDS(expTable_PROM_original,dim=c(3,4), col=batch)
####

# Boxplot
boxplot(select)
#boxplot(select[100:200])

# Media, desvio
means=colMeans(select)
std=apply(select, 2,  sd)
names=names(means)
c=1:length((names))

plot(c, means, type = "l", col="blue", xlab = "", ylab ="", ylim=range(c(-0.8, 0.8)) )
lines(c, std, col = "red")
lines(c, std*-1, col = "red")

plot(means, type = "l")
plot(std, type = "l")

 # Hetmap
# https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
df = t(select)
df <- scale(df)
heatmap(df, scale = "none")


