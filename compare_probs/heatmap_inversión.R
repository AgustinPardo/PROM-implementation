setwd("/Users/karinaatriztan/Documents/PTW-genome/Heatmaps") 

outpathcount2 = ("/Users/karinaatriztan/Documents/PTW-genome/Heatmaps/")


dir.create(outpathcount2, showWarnings=FALSE)

#BiocManager::install("ggplot2")
#BiocManager::install("ggdendro")
#BiocManager::install("reshape2")

#BiocManager::install("factoextra")
#BiocManager::install("purrr")
#BiocManager::install("dendextend")
library("ggplot2")
library("ggdendro")
library("reshape2")
library("grid")
library("factoextra")
library ("multiClust")
library("purrr")
library("dendextend")
#Load gene expression matrix
data.exprs <- read.table(file = "pathways_injury_fungivory_heatmap.txt", header = TRUE, sep="\t", comment.char="", row.names = 1)
str(data.exprs)
#promedios = t(apply(data.exprs, 1, function(col) c(mean(col[1:3]), mean(col[4:6]), mean(col[7:9]),  mean(col[10:10]), mean(col[11:12] , mean(col[13:14])))))
dev.off()
#colnames(data.exprs) = c("Alternaria_AC", "Alternaria_BC","Alternaria_C", "Plant_C","Plant_I3","Plant_I5")


#promedios = t(apply(data.exprs, 1, function(col) c(mean(col[1:5]), mean(col[6:8]), mean(col[9:11]),  mean(col[12:14]), mean(col[15:17]),  mean(col[18:20]),
                                                   #mean(col[21:23]) ,  mean(col[24:26]) ,  mean(col[27:29]))))
                                                  
#colnames(promedios) = c("Control", "MI30","MI90", "MI4h", "MI8h","Dm30","Dm90", "Dm4h", "Dm8h")
#Clustering de perfiles genéticos
#Unsupervised clustering and informative gene selection

#data.log <- log2(promedios)
#data.log[is.na(data.log)] <- 0


#This step  is for to scalate the data expresion because genes with high expression can affect the posterior analysis
data.exprs  <- as.matrix(data.exprs)
scaledata <- t(scale(t(data.exprs))) # Centers and scales data.
scaledata <- scaledata[complete.cases(scaledata),]

str(data.exprs)
# Once that data is scalated we performed the clusterization of rows and columns using the herarchical clustering + Pearson and Spearman correlation by cluster
hr <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete") # Clusters columns by Spearman correlation.

#te a heatmap of the data. 
#I like heatmap.2 in the gplots package, but you can pass the same call below to the base
#function heatmap() if you like:
library (RColorBrewer)
library(gplots)


gr.row <- cutree(hr, 6)

colr <- colorRampPalette(brewer.pal(6, "Set1"))(6)

tiff("MI_Dm_pathways_3.tiff", width = 8, height = 15, res= 100, units = "in")

colMain <- colorRampPalette(brewer.pal(9, "RdYlBu"))(100)
heatmap.2(scaledata,
          #Rowv= T,
          Rowv=as.dendrogram(hr), 
          #Colv= FALSE,
          Colv=as.dendrogram(hc),
          col=colMain,
          scale="none",
          margins = c(7, 15),
          cexCol = 1,
          key = T,
          #key.title = "GOterms",
          #key.xlab = NA,
          key.ylab = "gene number",
          keysize = 1,
          density.info = "none",
          labRow = rownames(scaledata),
          cexRow = 0.5,
          sepcolor = "gray",
          sepwidth = 0.002,
          colsep= 1:8, 
          #rowsep = 1:114,
          RowSideColors = colr[gr.row],
          trace = "none")

dev.off()


#cortar el dendrograma en clusters
library(pvclust)
hclusth1.5 = cutree(hr, h=1.5) #cut tree at height of 1.5
hclusth1.0 = cutree(hr, h=1.0) #cut tree at height of 1.0
hclusth0.5 = cutree(hr, h=0.5) #cut tree at height of 0.5

head (hclusth1.5)
str(hclusth1.0)

library(dendextend)
#plot the tree
plot(hr,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height",
     cex =0.1)
dev.off()
#add the three cluster vectors
the_bars <- cbind(hclusth0.5, hclusth1.0, hclusth1.5)
#this makes the bar
colored_bars(the_bars, hr, sort_by_labels_order = T, y_shift=-0.1, rowLabels = c("h=0.5","h=1.0","h=1.5"),cex.rowLabels=0.7)
#this will add lines showing the cut heights
abline(h=1.5, lty = 2, col="grey")
abline(h=1.0, lty = 2, col="grey")
abline(h=0.5, lty = 2, col="grey")

dev.off()

#Pretty neat! We see how cutting the tree at different heights
#yields different clusters. Alternatively, we can designate the number of clusters
#we want by using the ‘k’ option in cutree().
#Here I’m asking for four (k=4):
hclustk4 = cutree(hr, k=6)
plot(hr,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")
colored_bars(hclustk4, hr, sort_by_labels_order = T,
             y_shift=-0.1, rowLabels = c("k=6"),cex.rowLabels=2)
dev.off()
####################clustering OPTION 3################################ç
#Other option for obtain clusters, we can use the option cutreeDynamyc from the  "dynamicTreeCut " package. 
#Next method  cut the tree formed by hc using  the distM approach
library(dynamicTreeCut)
clusDyn <- cutreeDynamic(hr, distM = as.matrix(as.dist(1-cor(t(scaledata)))),
                         method = "hybrid")

plot(hr,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")

colored_bars(clusDyn, hr, sort_by_labels_order = T,
             y_shift=-0.1, rowLabels = NA,
             cex.rowLabels=0.7) #para 5660 genes produce 65 clusters diferentes
minModuleSize = 50;
#To this moment we opbserved that all clusering methos using produced around 65 clusterrs for out data


####################################################################


# Set the minimum module size
minModuleSize = 50;


###NOW WE GO TO  PERFORMED A CLUSTERIN USING THE WGCNA PACKAGE TO EXTRACT THE FORMED CLUSTERS TO FUTURE ANALYSIS
# Module identification using dynamic tree cut
#BiocManager::install("WGCNA")
library(WGCNA)

#WE WILL TO USE THE SAME DENDROGRAM BEFORE CONSTRUCTED AND WILL  CUTTED USING CUTREEDYNAMIC
dynamicMods = cutreeDynamicTree(dendro = hr, maxTreeHeight = 1.5, deepSplit = TRUE, minModuleSize = 25);
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
write.table(table(dynamicColors), file="keycolors4K.txt", sep="\t")
plotDendroAndColors(hr, dynamicColors, "k=6", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 1.5, main = "Gene dendrogram and module colors")


#discard the unassigned genes, and focus on the rest
restGenes= (dynamicColors != "grey")



#Extract modules # extraes el total de modules (cluster) que formaste
#y los guardas con respecto al color que  asisganaste, en nuestro paso, tenemos 65 K, los guarda sin el valor de expresión.
gene.names =row.names(data.exprs)
n=118
SubGeneNames =gene.names [1:n]
module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
  module= SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}


#########################################################################################
########################################################################################
#To calculate the cores (aka medoids) of each cluster we can use this function by 
#Biostar user Michael Dundrop. Here I’m using a cutheight of 1.5 to have a smaller 
#number of large clusters.

#Now we can plot the cores as a function of the samples.
#Since this is time series data, we use the time.point.hrs trait to plot the samples.
BiocManager::install("mudata")
library(mudata)
library(ggplot2)
library(reshape2)
#make a dataframe of cores by time

#Let’s perform the actual clsutering using K=4:

set.seed(20)
kClust <- kmeans(scaledata, centers=6, nstart = 1000, iter.max = 20)
kClusters <- kClust$cluster
write.table(kClust$cluster, "clustersize.txt", sep="\t", row.names=T, col.names=FALSE,quote=FALSE)

#Para ver el número de genes por cluster
kClust$centers

#Now we can calculate the cluster ‘cores’ aka centroids:

# function to find centroid in cluster i
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}
kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid, scaledata, kClusters)

#Plotting the centroids to see how they behave:
#library(ggplot2)
#BiocManager::install("reshape")
library(reshape)
#get in long form for plotting
Kmolten <- melt(kClustcentroids)
colnames(Kmolten) <- c('sample','cluster','value')

#plot
p1 <- ggplot(Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + 
  geom_line() +
  xlab("Fungivory vs injury") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Time",color = "Cluster")
p1


#So we have some interesting cluster profiles! If you do this analysis and recover
#cores that have very similar expression consider reducing your K.
#An a posteriori means of cluster validation is to correlate the cluster centroids with each other. 
#If the centroids are too similar then they will have a high correlation. 
#If your K number produces clusters with high correlation (say above 0.85) then consider reducing the number of clusters.
cor(kClustcentroids)
#To calculate the scores for a single cluster, in this case 2 we’ll extract the core data for cluster 2,
#then subset the scaled data by cluster =2. Then, we’ll calculate the ‘score’ by correlating each gene
#with the cluster core. We can then plot the results for each gene with the core overlayed:

#Subset the cores molten dataframe so we can plot the core
core2 <- Kmolten[Kmolten$cluster=="6",]

# #aca´guardar los genes pertenecientes a los clusters ,con su valor de expresion escaldo
K2 <- (scaledata[kClusters==6,])
str(K2)
write.table(K2, file="ID-cluster6_pathways_MIvsMd.txt", sep="\t")


#calculate the correlation with the core
corscore <- function(x){cor(x,core2$value)}
score <- apply(K2, 1, corscore)
#get the data frame into long format for plotting
K2molten <- melt(K2)#Acá podemos guardar la tabla
colnames(K2molten) <- c('gene','sample','value')
#add the score
K2molten <- merge(K2molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K2molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K2molten$order_factor <- 1:length(K2molten$gene)
#order the dataframe by score
K2molten <- K2molten[order(K2molten$gene),]
#set the order by setting the factors
K2molten$order_factor <- factor(K2molten$order_factor , levels = K2molten$order_factor)

# Everything on the same plot
p2 <- ggplot(K2molten, aes(x=sample,y=value), size = 4) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('gray','gray')) +
  theme_bw() +
  #this adds the core 
  geom_line(data=core2, aes(sample,value, group=cluster), color="orange",size =5, inherit.aes=FALSE) +
  xlab("Fungivory/Injury") +
  ylab("Expression") +
  labs(title= "Cluster 6",color = "Score", size = 15, color= "black") +
  theme(axis.text.x = element_text(size=rel(1)),
        axis.text.y.left =element_text(size=rel(3)),
        axis.title = element_text(size=rel(2.5)))
p2 

