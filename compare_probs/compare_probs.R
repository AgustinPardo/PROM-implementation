library (gplots)
library (RColorBrewer)

probtfgene_Sanz_437GEO <- read.csv("~/FBA_Tesis/PROM_trabajo/run_PROM_sets/probtfgene_Sanz_437GEO.txt")
probtfgene_Sanz_colombos_1021 <- read.csv("~/FBA_Tesis/PROM_trabajo/run_PROM_sets/probtfgene_Sanz_colombos_1021.txt")

probtfgene_Sanz_437GEO_r_t = paste(probtfgene_Sanz_437GEO$z_regulator,probtfgene_Sanz_437GEO$z_targets,sep="-")

# Combino las columnas de probabilidades en una matriz
mat <- cbind(probtfgene_Sanz_437GEO$probtfgene, probtfgene_Sanz_colombos_1021$probtfgene)
row.names(mat)= probtfgene_Sanz_437GEO_r_t
dim(mat)
colMain <- colorRampPalette(brewer.pal(9, "RdYlBu"))(100)
heatmap.2(mat, keysize=2, Colv = NA,  labCol = c("GEO437", "COLOMBOS"), cexCol = 1, trace="none", srtCol=0, key.title="", col=colMain )  

# dev.off()
# graphics.off()

# Remuevo filas si son las dos igual a 1
mat_no1=as.data.frame(mat)
mat_no1 <- mat_no1[!(mat_no1$V1 == 1 & mat_no1$V2 == 1 ),] 
dim(mat_no1)

# Hago histogramas
#Hao histograma con los que no son 1 en ambos, y saco ademas los 1 de cada set
mat_no1_GEO437=mat_no1[mat_no1$V1 !=1,]
mat_no1_GEO437=mat_no1_GEO437[,1]
mat_no1_COLOMBOS=mat_no1[mat_no1$V2 !=1,]
mat_no1_COLOMBOS=mat_no1_COLOMBOS[,2]

h = hist(mat_no1_GEO437, plot=FALSE, breaks=seq(0,1,by=0.05)) 
h$density = h$counts/sum(h$counts)
plot(h, xlab = "P", main= "Sanz GEO437",freq = FALSE, ylim=c(0,0.2),xlim=c(0,1))

h = hist(mat_no1_COLOMBOS, plot=FALSE, breaks=seq(0,1,by=0.05))
h$density = h$counts/sum(h$counts)
plot(h, xlab = "P", main= "Sanz COLOMBOS",freq = FALSE, ylim=c(0,0.2), xlim=c(0,1))

# Veo cuantos hay que tiene P=1 en cada set de datos
GEO437_1=mat_no1[mat_no1$V1 ==1,]
COLOMBOS_1=mat_no1[mat_no1$V2 ==1,]

mat_no1=as.matrix(mat_no1)
heatmap.2(mat_no1, keysize=2, Colv = NA,  labCol = c("GEO437", "COLOMBOS"), cexCol = 1, trace="none", srtCol=0, key.title="", col=colMain )  

