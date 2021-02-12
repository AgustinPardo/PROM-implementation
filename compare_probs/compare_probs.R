library (gplots)
library (RColorBrewer)

probtfgene_Sanz_437GEO <- read.csv("~/FBA_Tesis/PROM_trabajo/run_PROM_sets/probtfgene_Sanz_437GEO.txt")
probtfgene_Sanz_colombos_1021 <- read.csv("~/FBA_Tesis/PROM_trabajo/run_PROM_sets/probtfgene_Sanz_colombos_1021.txt")

probtfgene_Sanz_437GEO_r_t = paste(probtfgene_Sanz_437GEO$z_regulator,probtfgene_Sanz_437GEO$z_targets,sep="-")

# Combino las columnas de probabilidades en una matriz
mat <- cbind(probtfgene_Sanz_437GEO$probtfgene, probtfgene_Sanz_colombos_1021$probtfgene)
row.names(mat)= probtfgene_Sanz_437GEO_r_t

colMain <- colorRampPalette(brewer.pal(9, "RdYlBu"))(100)
heatmap.2(mat, keysize=2, Colv = NA,  labCol = c("GEO437", "COLOMBOS"), cexCol = 1, trace="none", srtCol=0, key.title="", col=colMain )  

# dev.off()
# graphics.off()

# Remuevo filas si son las dos igual a 1
mat_no1=as.data.frame(mat)
mat_no1 <- mat_no1[!(mat_no1$V1 == 1 & mat_no1$V2 == 1 ),] 

# Veo cuantos hay que tiene P=1 en cada set de datos
GEO437_1=mat_no1[!mat_no1$V1 !=1,]
COLOMBOS_1=mat_no1[!mat_no1$V2 !=1,]

mat_no1=as.matrix(mat_no1)
heatmap.2(mat_no1, keysize=2, Colv = NA,  labCol = c("GEO437", "COLOMBOS"), cexCol = 1, trace="none", srtCol=0, key.title="", col=colMain )  
