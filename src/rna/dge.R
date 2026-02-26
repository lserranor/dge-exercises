library("airway")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("EnhancedVolcano")
library("vsn")


## Cargamos nuestro dataset de ejemplo
data("airway")

summarized_experiment <- airway

## Vamos a explorar el colData, es decir, la info. experimental
exp_info <- as.data.frame(summarized_experiment@colData)

##Tambien miraremos que las cuentas de la matriz esten crudas, es decir que sean números enteros
matrix <- assay(summarized_experiment)

## Creamos el objeto DESeq2
dds <- DESeqDataSet(summarized_experiment, design = ~cell + dex)

## Prefiltrado. Vamos a eliminar genes
## con tan poquitas cuentas que no merece la pena conservar.
## Así ahorramos memoria, aunque no ganamos nada a efectos estadísticos.
keep <- rowSums(counts(dds)) >= 10 ## Seleccionar genes con más de 10 cuentas en todos los samples
dds <- dds[keep, ]

## Antes de hacer la DGE, hagamos un análisis exploratorio
## la función VST normaliza y estabiliza la varianza de las counts. Ideal para
## clustering o PCA.

## La diferencia entre VST y la función DESeq radica en que la expresión diferencial
## no usa las cuentas normalizadas (y estabilizadas) per sé, sino que los factores de 
## normalización y dispersión se incluyen en un modelo BN con las cuentas crudas.
## Para visualización, clustering etc. la función vst transforma las counts y estabiliza
## la varianza APARTE para que podamos trabajar con cuentas ya "listas".

vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "cell")
plotPCA(vsd, intgroup = "dex")

## Calculamos las distancias a partir de las cuentas normalizadas y variance-stabilized (vst)
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$cell, vsd$dex, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$cell, vsd$dex, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


## Nos preparamos para la DGE haciendo el modelo lineal a partir de la BN
## La función DESeq realiza todos los pasos de DESeq2 desde estimar los size factors
## hasta controla la dispersión
dds2 <- DESeq(dds, test = "Wald")

## Vamos a verificar cómo ha quedado la estimación de la dispersión 
plotDispEsts(dds2)

## Datos originales
msd <- meanSdPlot(counts(dds2), plot = FALSE)
msd$gg + ylim(0, 3000) ## Eliminamos los valores más altos sólo para la visualización

## Tras VST
meanSdPlot(assay(vsd))

## Exploremos cómo quedan los cambios de fold entre condiciones con respecto
## a las cuentas normalizadas
plotMA(dds2)

## Obtengamos nuestra lista de genes DEG
my_results <- results(object = dds2,
                      contrast = c("dex", "trt", "untrt"),
                      alpha = 0.05,
                      pAdjustMethod = "BH",
                      tidy = TRUE
                      )
#En la matriz de resultados nos aparecen varias columnas
# row: Nombre del gen en la anotación en la que se encontraban los genes en la matriz de expresión
# basseMean: media de cuentas normalizadas en todas las muestras
# log2FoldChange: La medida de cambio en log2. Positivo nos indicaría que está mas expresado en
# el primer termino de la comparativa (trt) y negativos más en el segundo término (untrt).
# Nota: EL logaritmo nos permite transformar los valores de Fold Change (trt/untrt) de una escala
# de 0 a Inf donde 1 es igualmente expresado en ambas condiciones >1 es mas expresado en trt
# y <1 mas expresado en un untrt a una escala de -Inf a -Inf donde 0 (log2(1) == 0) sería igualmente expresados
# en ambas condiciones. Para saber cuantas veces esta diferencialmente expresado debemos deshacer el log2
# FC = 2 ^ log2FC  
# lfcSE = Desviación estandar del log2FC
# stat= estadístico
# pvalue = p-valor simple
# padj = p-valor ajustado



##Anotamos los genes para mostrar los símbolos para facilitar el estudio biológico
genesID <- mygene::queryMany(my_results$row, scopes="ensembl.gene", fields="symbol", species="human")
genesID <- genesID[!duplicated(genesID$query),]
my_results$row <- ifelse(is.na(genesID$symbol),genesID$ query,genesID$symbol)

## Suele ser buena idea establecer un corte a priori de log fold.
## NOTA: En los resultados aparecen tambien genes con un log fold < al threshold
## El principal cambio es que ahora la prueba estadśitica no testea como hipotesis
# nula que el logFC sea != a 0 si no que testea que sea >= 1 o <= a -1

my_results_threshold <- results(object = dds2,
                                contrast = c("dex", "trt", "untrt"),
                                lfcThreshold = 1,
                                alpha = 0.05,
                                pAdjustMethod = "BH",
                                tidy = TRUE
                                )
##Anotamos los genes para mostrar los símbolos para facilitar el estudio biológico
genesID_threshold <-mygene::queryMany(my_results_threshold$row, scopes="ensembl.gene", fields="symbol", species="human")
genesID_threshold <- genesID_threshold[!duplicated(genesID_threshold$query),]
my_results_threshold$row <- ifelse(is.na(genesID_threshold$symbol),genesID_threshold$query,genesID_threshold$symbol)

## Heatmap de los genes TOP DGE por p-valor ajustado
mat <- assay(vsd)[head(order(my_results_threshold$padj), 30), ] 
pheatmap(mat)

# Creamos el Volcano plot 

EnhancedVolcano(my_results,
lab = my_results$row,
x = "log2FoldChange",
y = "padj",
title = "DEG treated vs untreated",
FCcutoff = 1,
pCutoff = 0.05,
subtitle = NULL,
boxedLabels = TRUE,
drawConnectors = TRUE,
labSize = 6.0)

