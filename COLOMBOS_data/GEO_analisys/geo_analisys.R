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

# exp_PROM_Original.txt
expTable_PROM_original <- read.table("exp_PROM_Original.txt", header=TRUE, row.names=1, sep=",")

# Boxplot
boxplot(expTable_PROM_original[1:100])
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
df <- scale(expTable_PROM_original[1:30])
heatmap(df, scale = "none")
