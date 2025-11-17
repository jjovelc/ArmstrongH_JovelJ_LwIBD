# https://rawgit.com/andreasmock/MetaboDiff/master/vignettes/MetaboDiff_tutorial.html
setwd('/Users/juanjovel/jj/data_analysis/heatherArmstrong/living_w_IBD/mixOmics/analyses_230913')

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
source("helper_functions.R")
library(MetaboDiff)

############################
# Functions
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Functions
############################

metabolites <- read.table('metabolomics_only_0-52w.tsv', header = T, row.names = 1, sep = '\t')


metabolites_t <- t(metabolites)

metabolites <- apply(metabolites, 2, RangeNorm)

metadata <- read.table('metadata_0-52w.tsv', header = T, row.names = 1, sep = '\t')
metadata$timePoint_severity <- paste(metadata$disease_severity, metadata$week, sep = '_') 
timePoint_severity_table <- data.frame(timePoint_severity = metadata$timePoint_severity)

logical_vector <- metadata$timePoint_severity %in% c("remission_0", "flare_52")
metadata_0remission_52flare <- metadata[logical_vector,] 
metabolites_0remission_52flare <- metabolites[logical_vector,]

met_t <- t(metabolites_0remission_52flare)
metabolites_ppm <- met_t / colSums(met_t) * 1000000


extMetaboLine <- function(df, key){
  metabo <- data.frame(metabolites_ppm[key, ])
  colnames(metabo) <- c("met_conc")
  return(metabo)
}

acetoacetate <- extMetaboLine(metabolites_ppm, "ACETOACETATE")
methionine   <- extMetaboLine(metabolites_ppm, "METHIONINE")
hippurate    <- extMetaboLine(metabolites_ppm, "HIPPURATE")
acetyl_methionine <- extMetaboLine(metabolites_ppm, "N_ACETYL_METHIONINE")
cytidine <- extMetaboLine(metabolites_ppm, "CYTIDINE2_3_CYCLICMONOPHOSPHATE")
guanidine <- extMetaboLine(metabolites_ppm, "GUANOSINE3_5_CYCLICMONOPHOSPHATE")
                              
makeBoxplot <- function(df, metadata){
  colnames(acetoacetate) <- c("met_conc")
  df$group <- metadata$timePoint_severity
  ggplot(df, aes(x=group, y=met_conc, fill=group)) +
    geom_boxplot()+
    scale_y_log10() +
    theme_classic()
}

makeBoxplot(acetoacetate, metadata_0remission_52flare)
makeBoxplot(methionine, metadata_0remission_52flare)
makeBoxplot(hippurate, metadata_0remission_52flare)
makeBoxplot(acetyl_methionine, metadata_0remission_52flare)
makeBoxplot(cytidine, metadata_0remission_52flare)
makeBoxplot(guanidine, metadata_0remission_52flare)


annotations <- read.table('metabolites_0remission_52flare_annotation.txt', 
                          header = T, row.names = 1, sep = '\t')

annotations <- read.table('b.txt', header = T, row.names = 1, sep = '\t')

writeFile <- function(df, fileName){
  write.table(df, fileName, sep = '\t', quote = F, row.names = T)
}



colData <- metadata_0remission_52flare
assay <- t(metabolites_0remission_52flare)
rowData <- annotations

writeFile(rowData, "rowdata.txt")
writeFile(colData, "metadata_0remission_52flare.txt")
writeFile(assay, "metabolites_0remission_52flare.txt")
writeFile(met_t, "metabolites_0remission_52flare-nonNorm.txt")

(met <- create_mae(assay,rowData,colData))

(met = knn_impute(met,cutoff=0.4))
(met <- remove_cluster(met,cluster=2))
(met <- normalize_met(met))

quality_plot(met,
             group_factor="timePoint_severity",
             label_colors=c("darkseagreen","dodgerblue"))


multiplot(
  pca_plot(met,
           group_factor="timePoint_severity",
           label_colors=c("tomato","dodgerblue")),
  tsne_plot(met,
            group_factor="timePoint_severity",
            label_colors=c("tomato","dodgerblue")),
  cols=2)

met = diff_test(met,
                group_factors = c("timePoint_severity","sex_f"))

str(metadata(met), max.level=2)

# Save original settings
op <- par(no.readonly=TRUE)

par(mfrow=c(1,2))
volcano_plot(met, 
             group_factor="timePoint_severity",
             label_colors=c("green","red"),
             dm_cutoff=0.5,
             p_adjust = FALSE)
volcano_plot(met, 
             group_factor="timePoint_severity",
             label_colors=c("green","red"),
             dm_cutoff=0.5,
             p_adjust = TRUE)

# Restore original settings
par(op)

par(mfrow=c(1,2))
volcano_plot(met, 
             group_factor="sex_f",
             label_colors=c("brown","orange"),
             p_adjust = FALSE)
volcano_plot(met, 
             group_factor="sex_f",
             label_colors=c("brown","orange"),
             p_adjust = TRUE)

Ttest_results <- met@metadata$ttest_timePoint_severity_remission_0_vs_flare_52

writeFile(Ttest_results, "Ttest_results_metabolomics_remission0-flare52.tsv")


met <- met %>%
  diss_matrix %>%
  identify_modules(min_module_size=5) %>%
  name_modules(pathway_annotation="SUB_CLASS") %>%
  calculate_MS(group_factors=c("timePoint_severity","sex_f"))


