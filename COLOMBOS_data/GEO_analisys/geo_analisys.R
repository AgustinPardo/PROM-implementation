# install.packages("BiocManager")
# install.packages("forcats")
# install.packages("stringr")
# install.packages("ggplot2")
# install.packages("ggrepel")
# install.packages("readr")
# install.packages("tidyr")
# install.packages("survminer")
# BiocManager::install("GEOquery")
# BiocManager::install("limma")
# BiocManager::install("pheatmap")
# BiocManager::install("org.Hs.eg.db")

#Recursos
#https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html#Filtering_lowly_expressed_genes
#https://www.stat.purdue.edu/bigtap/online/docs/Introduction_to_Microarray_Analysis_GSE15947.html
#https://sbc.shef.ac.uk/geo_tutorial/tutorial.nb.html
#https://rpubs.com/sharbie/462299

# exp_PROM_Original.txt
expTable_PROM_original <- read.table("exp_PROM_Original.txt", header=TRUE, row.names=1, sep=",")

# Boxplot
boxplot(expTable_PROM_original)
boxplot(expTable_PROM_original[100:200])

# Media, desvio
means=colMeans(expTable_PROM_original)
std=apply(expTable_PROM_original, 2,  sd)
names=names(means)
c=1:length((names))

plot(c, means, type = "l", col="blue", xlab = "", ylab ="", ylim=range(c(-0.8, 0.8)) )
lines(c, std, col = "red")
lines(c, std*-1, col = "red")

plot(means, type = "l")
plot(std, type = "l")

 # Hetmap
# https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
df = t(expTable_PROM_original)
df <- scale(df)
heatmap(df, scale = "none")
