setwd("/Users/karinaatriztan/Desktop/genomes/RNASEQ-PacBIo/HS02-2/DEG") 

outpathcount2 = "/Users/karinaatriztan/Desktop/genomes/RNASEQ-PacBIo/HS02-2/DEG/newHM-july2020/"


dir.create(outpathcount2, showWarnings=FALSE)

#BiocManager::install("ggplot2")
#BiocManager::install("ggdendro")
#BiocManager::install("reshape2")

library("ggplot2")
library("ggdendro")
library("reshape2")
library("grid")

otter <- read.csv(file = "all.txt", header = TRUE, sep="\t", comment.char="")
head (otter)
rownames(otter) <- NULL 

otter.scaled <- otter
otter.scaled[, c(2:8)] <- scale(otter.scaled[, 2:8])

 
 #Run clustering
otter.matrix <- as.matrix(otter.scaled)
rownames(otter.matrix) <- otter.scaled$gene.ID
otter.dendro <- as.dendrogram(hclust(d = dist(x = otter.matrix))) #Acá se corre el analisis de la distancia euclideana
 
 # Create dendrograma
dendro.plot <- ggdendrogram(data = otter.dendro, rotate = TRUE) +
        theme(axis.text.y = element_text(size = 0.2))
#tiff("dendro-all+1.tiff", width = 30, height = 40, res= 100, units = "in") #Para exportar la imagen en alta calidad
 # Preview the plot
print(dendro.plot)
dev.off() 

# Data wrangling
otter.long <- melt(otter.scaled, id = c("gene.ID"))
# Extract the order of the tips in the dendrogram
otter.order <- order.dendrogram(otter.dendro)
# Order the levels according to their position in the cluster
otter.long$gene.ID <- factor(x = otter.long$gene.ID,
                               levels = otter.scaled$gene.ID[otter.order], 
                               ordered = TRUE)
# Create heatmap plot
tiff("HM3-3.tiff", width = 25, height = 35, res= 100, units = "in") #Para exportar la imagen en alta calidad

heatmap.plot1 <- ggplot(data = otter.long, aes(x = variable, y = gene.ID, fill= value )) +
        geom_tile()+
        scale_fill_gradient2(low = ("brown4"),#chartreuse
                             mid = "white",
                             high = ("blue4"),
                             limits = c( -5,5), 
                             breaks = c(-2, -1, 3, 7),
                             midpoint = 0,
                             ) +
        theme(axis.text.y = element_text(size = 0),
              axis.text.x = element_text(size = 50),
              legend.position = "left", legend.key.size = unit(5, "cm"),
              legend.text = element_text(size = 30),
              legend.title = element_text(size =40)
              )
              

print(heatmap.plot1)

dev.off() 

#Para exportar la imagen en alta calidad

#tiff("HM+dendro.TIFF", width = 30, height = 40, res= 100, units = "in")
#grid.newpage()
#print(heatmap.plot1, vp = viewport(x = 4, y = 0.5, width = 0.8, height = 0.9))
#print(dendro.plot, vp = viewport(x = 4, y = 0.43, width = 0.2, height = 1.0))

dev.off() 


#Para extracción de clusters formados
otter.dendro <- as.dendrogram(hclust(d = dist(x = otter.matrix)))

# Create dendro
dendro.plot <- ggdendrogram(data = otter.dendro, rotate = TRUE) +
        theme(axis.text.y = element_text(size = 0.2))
#  
library (ggplot2)



otter.matrix <- as.matrix(otter.scaled)
rownames(otter.matrix) <- otter.scaled$gene.ID
otter.dendro <- as.dendrogram(hclust(d = dist(x = otter.matrix)))

# Create dendro
dendro.plot <- ggdendrogram(data = otter.dendro, rotate = TRUE) +
        theme(axis.text.y = element_text(size = 0.2))


#Make a tSNE plot Note: tSNE is a stochastic method. Everytime you run it you will get slightly different results. 
#For convenience we can get the same results if we seet the seed.
#BiocManager::install("Seurat")

#Para extraer los cluster, usaremos cutree
distance_eucledian <- dist(otter.matrix, method= "euclidean")

#Perform hierarchical clustering using ward linkage
summary (otter)
plot (I30 ~ gene.ID ,  otter, cex =0.2)
