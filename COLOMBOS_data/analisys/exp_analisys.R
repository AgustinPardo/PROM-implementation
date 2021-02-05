install.packages("BiocManager")
install.packages("forcats")
install.packages("stringr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readr")
install.packages("tidyr")
install.packages("survminer")
BiocManager::install("GEOquery")
BiocManager::install("pheatmap")
BiocManager::install("org.Hs.eg.db")

#limma
BiocManager::install("limma")
#quantile
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

#colombos_GPL <- read.table("/home/agustin/FBA_Tesis/PROM_trabajo/COLOMBOS_data/colombos_mtube_exprdata_20151029.txt",  dec = "." , skip=5,  header=FALSE, fill=TRUE)
colombos_mtube_GPL <- read.delim("~/FBA_Tesis/PROM_trabajo/COLOMBOS_data/colombos_mtube_GPL.txt", header=FALSE)
colombos_mtube_GPL <- as.matrix (colombos_mtube_GPL)
#colombos_GPL <- colombos_GPL[1,]
batch=colombos_mtube_GPL
batch=as.vector(batch)
a<-c(1,2,3)
colombos_all <- read.table("/home/agustin/FBA_Tesis/PROM_trabajo/COLOMBOS_data/colombos_all.txt", row.names=1, dec = "." ,skip=7,  header=FALSE, fill=TRUE)
select <- as.matrix (colombos_all)
boxplot(select,outline=FALSE)

#select[is.na(select)] <- 0
#select <- as.matrix (colombos_9776)

# limma 
library("limma")
select_quantile <- normalizeBetweenArrays(select, method="quantile")
boxplot(select_quantile,outline=FALSE)

#quantile
#library("preprocessCore")
#select_quantile=normalize.quantiles(select)
#boxplot(select_quantile,outline=FALSE)

# # rma
# library("metafor")
# data.rma = rma(colombos_GSM27855)
# data.matrix = exprs(data.rma)
# 
# # normalize
# library("EnrichmentBrowser")
# maSE <- normalize(colombos_GSM27855) 

#ComBat

# a=replicate(138, 1)
# b=replicate(299, 2)
# # a=replicate(35, 1)
# # b=replicate(2, 2)
# # c=replicate(17, 3)
# # d=replicate(16, 4)
# # e=replicate(6,5)
# batch<- c(a,b)
# length(batch)
# dim(select_quantile)

library("sva")
select_quantile[is.na(select_quantile)] = 0
data_combat=ComBat(dat=select_quantile, batch=batch, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
boxplot(data_combat,outline=FALSE)


#####
# MDS con combat
mds = plotMDS(data_combat,xlim=c(-2,2) ,xlab = "Dim1", ylim=c(-2,2), ylab = "Dim2", cex.lab= 1.5,pch = 15,cex=2, col=c( rep("black",6),  rep("red",12), rep("yellow",12), rep("orange",12))) #,rep("green",3),rep("deeppink",3),rep("gray",3), rep("pink",3)))
batch_col=as.numeric(factor(batch))
plotMDS(data_combat,dim=c(3,4), col=batch_col, ylim=c(-0.6, 0.6), xlim=c(-0.6, 0.6))
xx = c(-0.4,-0.4,0.4,0.4)
yy = c(-0.4,0.4,0.4,-0.4)
# When density=0, col refers to the line colour
polygon(xx,yy, density=0, col="black")

# MDS sin combat
mds = plotMDS(select_quantile,xlim=c(-2,2) ,xlab = "Dim1", ylim=c(-2,2), ylab = "Dim2", cex.lab= 1.5,pch = 15,cex=2, col=c( rep("black",6),  rep("red",12), rep("yellow",12), rep("orange",12))) #,rep("green",3),rep("deeppink",3),rep("gray",3), rep("pink",3)))
batch_col=as.numeric(factor(batch))
plotMDS(select_quantile, dim=c(3,4), col=batch_col ,ylim=c(-0.6, 0.6), xlim=c(-0.6, 0.6))
xx = c(-0.4,-0.4,0.4,0.4)
yy = c(-0.4,0.4,0.4,-0.4)
# When density=0, col refers to the line colour
polygon(xx,yy, density=0, col="black")

#Normalizamos de nuevo
select_quantile2 <- normalizeBetweenArrays(data_combat, method="quantile")
boxplot(select_quantile2,outline=FALSE)

# Set final. Exporto la matriz


####
# 
# # Boxplot
# boxplot(select)
# #boxplot(select[100:200])
# 
# select=select_quantile
# select[is.na(select)] = 0
# # Media, desvio
# means=colMeans(select)
# std=apply(select, 2,  sd)
# names=names(means)
# c=1:length((names))
# 
# plot(c, means, type = "l", col="blue", xlab = "", ylab ="", ylim=range(c(-0.8, 0.8)) )
# lines(c, std, col = "red")
# lines(c, std*-1, col = "red")
# 
# plot(means, type = "l")
# plot(std, type = "l")
# 
#  # Hetmap
# # https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
# df = t(select)
# df <- scale(df)
# heatmap(df, scale = "none")