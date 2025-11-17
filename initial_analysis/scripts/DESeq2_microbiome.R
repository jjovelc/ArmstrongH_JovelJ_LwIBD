###########################################################################
#         This script conducts differential expression analysis           #
#         using the DESeq2 R package. In addition, it generates           #
#         exploratory plots (hierarchical clustering and PCoA).           #
#         It exports results with significant results to an Excel         # 
#         file. The same set of data, plus normalized input data          #
#         is exported to a tsv file. Finally, it depicts DE               #
#         features in a volcano plot.                                     #
#                                                                         #
#         Author: Juan Jovel (juan@dayhoff.ai)                            #
###########################################################################

# To install DESeq2, run the following commands:
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

# All other required libraries:
# "xlsx", "ecodist", "ggplot2"... etc
# can be installed with command:
# install.packages("library") [for individual libraries]
# or: install.packages(c("library1","library2","libraryn")) [for multiple libraries]

library(DESeq2)
library(xlsx) 
library(ecodist)
library(RColorBrewer)
library(pheatmap)
library(vegan)
library(ggplot2)

setwd("/Users/juanjovel/jj/data_analysis/heatherArmstrong/living_w_IBD/mixOmics/analyses_230913")
data_file     <- "taxa_only_0-52w_RF_counts.tsv"
metadata_file <- "taxa_only_0-52w_RF_meta.tsv"


heatmap_file <- gsub('.tsv', '_HCheatmap.pdf', data_file)
PCoA_file    <- gsub('.tsv', '_PCoA.pdf', data_file)
volcano_file <- gsub('.tsv', '_volcanoPlot.pdf', data_file)

allRes_file       <- gsub('.tsv', '_allRes_normData.tsv', data_file)
resSign_xlsx_file <- gsub('.tsv', '_sign_DESeq2_res.xlsx', data_file)
res0.1_xlsx_tab   <- gsub('.tsv', '_p0.1_DE_features', data_file)
res0.05_xlsx_tab  <- gsub('.tsv', '_p0.05_DE_features', data_file)
res0.01_xlsx_tab  <- gsub('.tsv', '_p0.01_DE_features', data_file)

# Read data
data <- read.table(data_file, header = T, row.names = 1, sep = "\t")


metadata <- read.table(metadata_file, header = T, row.names = 1, sep = "\t")

# Remove decimals because DESeq2 only handles integer
data <- round(as.matrix(data, digits = 0))

# Import data into DESEq2 object and define your design
ddsMat <- DESeqDataSetFromMatrix(countData = data, colData = metadata, 
                                 design =~ group)

# Apply logarithmic transformation to create exploratory plots 
# (Hierarchical clustering and PCoA)
rld <- rlog(ddsMat, fitType = 'local')

# Hierarchical clustering
bc_dist <- bcdist(t(assay(rld)))
sampleDistMatrix <- as.matrix(bc_dist)
rownames(sampleDistMatrix) <- paste(rld$label, sep = '_')
colnames(sampleDistMatrix) <- paste(rld$label, sep = '_')

colors <- colorRampPalette(rev(brewer.pal(9, 'RdBu')))(255)

# Generate heatmap
pdf(file=heatmap_file, width = 10, height = 8)
heatmap <- pheatmap(sampleDistMatrix, color = colors, 
                    clustering_distance_cols = bc_dist,
                    clustering_distance_rows = bc_dist)
dev.off()

# Generate Principal Coordinates (PCoA)
PCo <- pco(bc_dist)

# Put principal coordinates into a dataframe

PCo_df <- data.frame(PCo$vectors[,1],PCo$vectors[,2])
col_names <- c("PCoA1", "PCoA2")
colnames(PCo_df) <- col_names

# Plot PCoA
pdf(file=PCoA_file, width = 10, height = 8)
p <- ggplot(PCo_df, aes(PCoA1, PCoA2, color=metadata$group)) +
  geom_point(shape=19, size=7) +
  scale_color_manual(values=c("firebrick1", "dodgerblue")) + 
  geom_text(aes(label=row.names(PCo_df)), color='black') +
  theme_bw()
print(p)
dev.off()

# Transform data using the harmonic mean
gm_mean <- function(x, na.rm=T) {
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMean <- apply(counts(ddsMat), 1, gm_mean)
ddsMat <- estimateSizeFactors(ddsMat, geoMean = geoMean)

ddsMat  <- DESeq(ddsMat, fitType = 'local')
res    <- results(ddsMat)
resultsNames(ddsMat)

resP0.1  <- subset(res, padj < 0.1)
resP0.05 <- subset(res, padj < 0.05)
resP0.01 <- subset(res, padj < 0.01)

# Bind row names (taxa) to results
resP0.1  <- cbind(Taxa=rownames(resP0.1), resP0.1)
resP0.05 <- cbind(Taxa=rownames(resP0.05), resP0.05)
resP0.01 <- cbind(Taxa=rownames(resP0.01), resP0.01)

write.xlsx(resP0.1, file=resSign_xlsx_file, sheetName = res0.1_xlsx_tab, append=FALSE)
write.xlsx(resP0.05, resSign_xlsx_file, sheetName = res0.05_xlsx_tab, append=TRUE)
write.xlsx(resP0.01, resSign_xlsx_file, sheetName = res0.01_xlsx_tab, append=TRUE)

write.table(resP0.1, paste(res0.1_xlsx_tab, 'tsv', sep='.'), sep = "\t", row.names = F, quote = F)
write.table(resP0.05, paste(res0.05_xlsx_tab, 'tsv', sep='.'), sep = "\t", row.names = F, quote = F)
write.table(resP0.01, paste(res0.01_xlsx_tab, 'tsv', sep='.'), sep = "\t", row.names = F, quote = F)

# Include normalized data in results table

resdata <- merge(as.data.frame(res), as.data.frame(counts(ddsMat, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Feature"
resdata <- cbind(geneName=rownames(resdata), resdata)
write.table(resdata, file=allRes_file, sep="\t", quote = F, row.names = F)


## Volcano plot with "significant" genes colored by significance
volcanoplot <- function(res, lfcthresh=2, sigthresh=0.1, main='', legendpos="topright", labelsig=TRUE, textcx=1) {
  max_lfc <- max(abs(res$log2FoldChange))
  max_padj <- max(-log10(res$padj))
  
  with(res, plot(log2FoldChange, -log10(padj), pch=20, main=main, cex=0.8, ylim=c(0, max_padj), xlim=c(-max_lfc, max_lfc)))
  with(subset(res, padj<sigthresh), points(log2FoldChange, -log10(padj), pch=20, col="dodgerblue", cex = 1.2))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(padj), pch=20, col="darkgrey", cex = 1.2))
  with(subset(res, padj<sigthresh & log2FoldChange > lfcthresh), points(log2FoldChange, -log10(padj), pch=20, col="red", cex = 1.2))
  with(subset(res, padj<sigthresh & log2FoldChange < -lfcthresh), points(log2FoldChange, -log10(padj), pch=20, col="forestgreen", cex = 1.2))
}

pdf(file=volcano_file, width=10, height=8)
volcanoplot(res, lfcthresh=1, sigthresh=0.05, textcx=0.5)

mycol <- c("black", "dodgerblue", "darkgray", "red", "forestgreen")
legend("topleft", legend=c("Padj >= 0.05", "Padj < 0.05; FC < 2", "Padj > 0.0.5; FC > 2", "Padj < 0.05; FC > 2", "Padj < 0.05; FC < -2"),
       pch=c(19, 19), col=mycol, bg=mycol, pt.cex=1, bty="n", cex=1, text.col="black", horiz=FALSE, inset=c(0.02, 0.02))
dev.off()
