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


# Broad GSEA uses human gene sets and human gene names even for gene sets generated from mouse/rat data 
# To use it for different species, gene names have to be converted to human orthologs.
# Note: download Mouse gene name remapping table from Broad website. Current verison: MsigDB v.7.0.
mouse2human <- read.table(file="./GSEA/Mouse_Gene_Symbol_Remapping_MSigDB.v7.0.chip", sep="\t", header=T, fill=T, quote="") %>% .[,1:2] %>% setNames(c("mouse",'human'))

output_sol <- merge(as.data.frame(sol), mouse2human, by.x='row.names', by.y='mouse', sort = F)
output_edl <- merge(as.data.frame(edl), mouse2human, by.x='row.names', by.y='mouse', sort = F)

colnames(output_sol)[1] <- 'Mouse_gene_symbol'
colnames(output_sol)[8] <- 'Human_ortholog'
colnames(output_edl)[1] <- 'Mouse_gene_symbol'
colnames(output_edl)[8] <- 'Human_ortholog'

for(i in 1:nrow(output_sol)) {
  pval <- output_sol$pvalue[i]
  if(pval %in% output_sol$pvalue[duplicated(output_sol$pvalue)]) { 
    iteration <- 0
    repeat {
      iteration <- iteration + 1
      pval <- pval + iteration * (output_sol$pvalue[i]/10**14)
      if (!(pval %in% output_sol$pvalue)) {
        output_sol$pvalue[i] <- pval
        break
      }
    }
  }
}

for(i in 1:nrow(output_edl)) {
  pval <- output_edl$pvalue[i]
  if(pval %in% output_edl$pvalue[duplicated(output_edl$pvalue)]) { 
    iteration <- 0
    repeat {
      iteration <- iteration + 1
      pval <- pval + iteration * (output_edl$pvalue[i]/10**14)
      if (!(pval %in% output_edl$pvalue)) {
        output_edl$pvalue[i] <- pval
        break
      }
    }
  }
}

# write rank files for GSEA
output_sol$rank <- -sign(output_sol$log2FoldChange) * log10(output_sol$pvalue)
output_edl$rank <- -sign(output_edl$log2FoldChange) * log10(output_edl$pvalue)

rank_sol <- output_sol[, c('Human_ortholog', 'rank')] %>% setNames(c('NAME','RANK'))
rank_edl <- output_edl[, c('Human_ortholog', 'rank')] %>% setNames(c('NAME','RANK'))
rank_sol <- rank_sol[order(rank_sol$RANK),]
rank_edl <- rank_edl[order(rank_edl$RANK),]
#sanity check# length(unique(rank$RANK)) == length(rank$RANK)
write.table(rank_sol, file="./GSEA/sol_gsea.rnk", row.names = F, quote = F, sep = "\t")
write.table(rank_edl, file="./GSEA/edl_gsea.rnk", row.names = F, quote = F, sep = "\t")
















