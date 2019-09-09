library(magrittr)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(dplyr)
library(stringr)
library(pheatmap)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))


#------------------------------------------- Import data ----------------------------------------------------
sampleTable <- read.table(file = "sampleTable.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
counts <- read.table(file = "geneCounts.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)


#-------------------------------------- PCA plots --------------------------------------------
# all samples
PlotData <- counts[,2:ncol(counts)]
PlotData <- PlotData[rowMeans(PlotData) >= 10, ] %>% as.matrix() %>% rlog()
PlotData <- PlotData[order(rowVars(PlotData), decreasing = TRUE), ] %>% .[1:1000,]
pca <- prcomp(t(PlotData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAsummary <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                         group = paste0(sampleTable$age, "-", sampleTable$tissue),
                         name  = sapply(names(pca$x[, 1]), function(x) { return(substr(x,2,100))   })
)

png(file = "PCA_all.png", width = 2400, height = 2400, units = 'px', res = 300)
ggplot(data = PCAsummary, aes_string(x = "PC1", y = "PC2", color = "group")) + 
  geom_point(size = 6) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(linetype = "solid", fill = NA, size = 1.5)) +
  theme(axis.ticks.length = unit(2, "mm"), axis.ticks = element_line(size = 1), axis.text = element_text(size = rel(1.25))) +
  theme(axis.title = element_text(size = rel(1.5)))+
  theme(legend.text = element_text(size = rel(1.2)), legend.title = element_blank())
dev.off()


# samples exluding a single major outlier
PlotData <- counts[,2:ncol(counts)]
PlotData <- PlotData[, colnames(PlotData) != 'X32' ]
PlotData <- PlotData[rowMeans(PlotData) >= 10, ] %>% as.matrix() %>% rlog()
PlotData <- PlotData[order(rowVars(PlotData), decreasing = TRUE), ] %>% .[1:1000,]
pca <- prcomp(t(PlotData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAsummary <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                         group = paste0(sampleTable[sampleTable$sampleID != '32', 'age'], "-", sampleTable[sampleTable$sampleID != '32', 'tissue']),
                         name  = sapply(names(pca$x[, 1]), function(x) { return(substr(x,2,100))   })
)

png(file = "test.png", width = 2400, height = 2400, units = 'px', res = 300)
ggplot(data = PCAsummary, aes_string(x = "PC1", y = "PC2", color = "group")) + 
  geom_point(size = 6) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(linetype = "solid", fill = NA, size = 1.5)) +
  theme(axis.ticks.length = unit(2, "mm"), axis.ticks = element_line(size = 1), axis.text = element_text(size = rel(1.25))) +
  theme(axis.title = element_text(size = rel(1.5)))+
  theme(legend.text = element_text(size = rel(1.2)), legend.title = element_blank())
dev.off()


###############################################################################################################################
###################################   Differential Gene Expression Analysis mRNAseq   #########################################
###############################################################################################################################
df <- counts[, colnames(counts) != 'X32' ]
row.names(df) <- df$geneName
df <- df[,2:ncol(df)]

colData <- data.frame(sample = colnames(df),
                      age = factor(sampleTable[sampleTable$sampleID != '32', 'age'], levels=c('young','old')),
                      tissue = factor(sampleTable[sampleTable$sampleID != '32', 'tissue'], levels=c('sol','edl')),
                      stringsAsFactors = FALSE
)

dds <- DESeqDataSetFromMatrix(countData = df, colData = colData, design =~ tissue + age + tissue:age)
dds <- dds[rowMeans(counts(dds)) >= 10, ]
dds <- DESeq(dds)
resultsNames(dds)


sol <- results(dds, name = "age_old_vs_young", cooksCutoff = FALSE, independentFiltering = FALSE)
sol <- sol[order(sol$padj),]

edl <- results(dds, contrast = list(c("age_old_vs_young","tissueedl.ageold" )), cooksCutoff = FALSE, independentFiltering = FALSE)
edl <- edl[order(edl$padj),]

ageDiff_edl_vs_sol <- results(dds, name = "tissueedl.ageold", cooksCutoff = FALSE, independentFiltering = FALSE)
ageDiff_edl_vs_sol <- ageDiff_edl_vs_sol[order(ageDiff_edl_vs_sol$padj),]

write.csv(sol, file="SOL_old_vs_young.csv")
write.csv(edl, file="EDL_old_vs_young.csv")
write.csv(ageDiff_edl_vs_sol, file="ageDiff_edl_vs_sol.csv")

















