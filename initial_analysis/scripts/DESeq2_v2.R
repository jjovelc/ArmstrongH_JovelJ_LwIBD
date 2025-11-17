#This script performs expolatory analysis and differential expression analysis
#on shotgun metagenomics datasets.
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(EnhancedVolcano)
library(DataEditR)
library(pheatmap)
library(ecodist)

setwd('/Users/juanjovel/jj/data_analysis/heatherArmstrong/living_w_IBD')

data_file <- 'all_samples_mpa_tax_shortNames.tsv'
metadata_file <- 'metadata.tsv'

##### Functions #####
makeHCheatmap <- function(rld, metadata, intgroup){
  # Bray-Curtis distances
  dist        <- bcdist(t(assay(rld)))
  dist_matrix <- as.matrix(dist)
  row.names(dist_matrix) <- metadata[[intgroup]]
  df_name <- deparse(substitute(metadata))
  
  # Use gsub to create the new name
  file_name <- "HC_plot.pdf"
  pdf(file=file_name, width = 10, height = 8)
  pheatmap(dist_matrix, 
           color=colorRampPalette(brewer.pal(n=9,name = "RdYlBu"))(255), 
           clustering_distance_cols = dist, 
           clustering_distance_rows = dist
  )
  dev.off()
}


# Make PCoA plot
makePCoA <- function(rld, metadata, intgroup){
  file_name <- "PCoA_plot.pdf"
  bc_dist <- bcdist(t(assay(rld)))
  PCo <- pco(bc_dist)
  PCo_df <- data.frame(PCo$vectors[,1],PCo$vectors[,2])
  col_names <- c("PCoA1", "PCoA2")
  colnames(PCo_df) <- col_names
  
  pdf(file=file_name, width = 10, height = 8)
  p <- ggplot(PCo_df, aes(PCoA1, PCoA2, color=metadata[[intgroup]])) +
    geom_point(shape=19, size=7) +
    geom_text(aes(label=row.names(PCo_df)), color='black') +
    theme_bw()
  print(p)
  dev.off()
}
##### Functions END #####

data <-read.table (data_file, header=T, sep = '\t', row.names = 1)
metadata <-read.table (metadata_file, header=T, sep = '\t', row.names = 1)


# Remove sample X299_wk52, which seems to be an ourlier
metadata_noOutlier <- subset(metadata, row.names(metadata) != "X299_WK52" )
data_noOutlier     <- data[, colnames(data) != "X299_WK52"]


#############
# Extract subset wk0 flare and wk52 remission
# 2. Create the logical vectors

wk0_flare <- metadata_noOutlier$timePoint == "wk0"  & metadata_noOutlier$severity == "flare"
wk52_remission <- metadata_noOutlier$timePoint == "wk52" & metadata_noOutlier$severity == "remission"

wk0_flare_meta <- subset(metadata_noOutlier, metadata_noOutlier$timePoint == "wk0" & metadata_noOutlier$severity == "flare")
wk52_remission_meta <- subset(metadata_noOutlier, metadata_noOutlier$timePoint == "wk52" & metadata_noOutlier$severity == "remission")
wk0_flare_wk52_remission_meta <- rbind(wk0_flare_meta, wk52_remission_meta)

# 3. Extract sample names based on the logical vectors
wk0_flare_samples <- data_noOutlier[,wk0_flare]
wk52_remission_samples <- data_noOutlier[,wk52_remission]
wk0_flare_wk52_remission_data <- cbind(wk0_flare_samples, wk52_remission_samples)

# Combine timePoint and severity in a single variable
wk0_flare_wk52_remission_meta$combinedStatusTime <- with(wk0_flare_wk52_remission_meta,
                                                         interaction(severity, timePoint, sep="_"))

ddsMatPat <- DESeqDataSetFromMatrix(countData = wk0_flare_wk52_remission_data, 
                                    colData = wk0_flare_wk52_remission_meta, 
                                    design =~ patient + combinedStatusTime)

keep      <- rowSums(counts(ddsMatPat) >= 10) >= 10
ddsMatPat <- ddsMatPat[keep,]

rldPat <- rlog(ddsMatPat, blind=FALSE)
makeHCheatmap(rldPat, wk0_flare_wk52_remission_meta, 'severity')
makePCoA(rldPat, wk0_flare_wk52_remission_meta, 'severity')
###############


# Combine timePoint and severity in a single variable
metadata$combinedStatusTime <- with(metadata, interaction(severity, timePoint, sep="_"))

ddsMatPat <- DESeqDataSetFromMatrix(countData = data, 
                                 colData = metadata, 
                                 design =~ patient + combinedStatusTime)

keep      <- rowSums(counts(ddsMatPat) >= 100) >= 10
ddsMatPat <- ddsMatPat[keep,]

                                 
# Transform data using the harmonic mean
gm_mean <- function(x, na.rm=T) {
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMean <- apply(counts(ddsMatPat), 1, gm_mean)
ddsMatPat <- estimateSizeFactors(ddsMatPat, geoMean = geoMean)


rldPat <- rlog(ddsMatPat, blind=FALSE)

makeHCheatmap(rldPat, metadata, 'severity')
makePCoA(rldPat, metadata, 'diagnosis')


#####################################
#     Include Patient's factor      #
#####################################


# Running DESeq with the LRT:
ddsMatPat <- DESeq(ddsMatPat, test="LRT", reduced=~ patient, fitType='local')
ddsMatPat2 <- nbinomLRT(ddsMatPat, full = design(ddsMatPat), 
                        type = c("DESeq2", "glmGamPoi"),
                        reduced=~ patient, maxit = 2000)

# Inspect the non-converging features
which(!mcols(ddsMatPat)$fullBetaConv)


res1 <- results(ddsMatPat, contrast=c("combinedStatusTime", "flare_wk0", "remission_wk52"), test="Wald")
sign_res1 <- subset(res1, padj < 0.05)
nrow(sign_res1)
write.table(sign_res1, "flare-wk0_vs_remission-wk52_pAdj0.05.tsv",
            quote = F, sep = '\t')


res2 <- results(ddsMatPat, contrast=c("combinedStatusTime", "remission_wk0", "flare_wk52"), test="Wald")
sign_res2 <- subset(res2, padj < 0.05)
nrow(sign_res2)
write.table(sign_res2, "remission-wk0_vs_flare-wk52_pAdj0.05.tsv",
            quote = F, sep = '\t')


res4 <- results(ddsMatPat, contrast=c("combinedStatusTime", "flare_wk0", "remission_wk0"), test="Wald")
sign_res4 <- subset(res4, padj < 0.05)
nrow(sign_res4)

res5 <- results(ddsMatPat, contrast=c("combinedStatusTime", "flare_wk26", "remission_wk26"), test="Wald")
sign_res5 <- subset(res5, padj < 0.05)
nrow(sign_res5)

res6 <- results(ddsMatPat, contrast=c("combinedStatusTime", "flare_wk52", "remission_wk52"), test="Wald")
sign_res6 <- subset(res6, padj < 0.05)
nrow(sign_res6)


#####################################
#      Create sub-dataframes        #
#####################################
# Modify the format of WK26 samples names
new_metadata <- metadata %>%
  mutate(sample_name = if_else(
    timePoint == "wk26",
      as.numeric(str_extract(row.names(metadata), "(?<=W26_P)\\d+")),
    as.numeric(str_extract(row.names(metadata), "\\d+"))
  ))


# Extracting Patients that were in remission at wk0 
# and turned to flare at wk26 (RF_wk0_26)

# Subset 1
#################################
joined_data <- new_metadata %>%
  filter(timePoint == "wk0" & severity == "remission") %>%
  inner_join(new_metadata %>% filter(timePoint == "wk26" & severity == "flare"), by = "sample_name")

# Splitting data for wk0 and wk26 and formatting the sample names
wk0_data <- joined_data %>%
  transmute(sample_name = paste0("X", sample_name, "_Baseline"),
            severity = severity.x)

wk26_data <- joined_data %>%
  transmute(sample_name = paste0("W26_P", sample_name),
            severity = severity.y)

# Binding rows
RF_wk0_26 <- bind_rows(wk0_data, wk26_data)
row.names(RF_wk0_26) <- RF_wk0_26$sample_name

# Extracting the desired columns from counts table
selected_counts <- data %>%
  select(all_of(RF_wk0_26$sample_name))

# Print the result
RF_wk0_26$sample_name <- NULL

data_matrix <- as.matrix(round(selected_counts, digits=0))
# Import data into DESEq2 object and define your design
ddsMat <- DESeqDataSetFromMatrix(countData = data_matrix, colData = RF_wk0_26, 
                                 design=~severity)

keep <- rowSums(counts(ddsMat) >= 10) >= 3
ddsMat <- ddsMat[keep,]

# Transform data using the harmonic mean
geoMean <- apply(counts(ddsMat), 1, gm_mean)
ddsMat  <- estimateSizeFactors(ddsMat, geoMean = geoMean)

rld <- rlog(ddsMat, blind=FALSE, fitType = 'local')

makeHCheatmap(rld, RF_wk0_26, severity)
makePCoA(rld, RF_wk0_26, severity)

ddsMat <- DESeq(ddsMat, fitType='local')
resultsNames(ddsMat)
res <- results(ddsMat, contrast = c("severity", "flare", "remission"))
sign_res <- subset(res, padj < 0.05)
nrow(sign_res)

# Subset 2
##################################
# Extracting Patients that were in flare at wk0 
# and turned to remission at wk26 (FR_wk0_26)
joined_data <- new_metadata %>%
  filter(timePoint == "wk0" & severity == "flare") %>%
  inner_join(new_metadata %>% filter(timePoint == "wk26" & severity == "remission"), by = "sample_name")

# Splitting data for wk0 and wk26 and formatting the sample names
wk0_data <- joined_data %>%
  transmute(sample_name = paste0("X", sample_name, "_Baseline"),
            severity = severity.x)

wk26_data <- joined_data %>%
  transmute(sample_name = paste0("W26_P", sample_name),
            severity = severity.y)

# Binding rows
FR_wk0_26 <- bind_rows(wk0_data, wk26_data)
row.names(FR_wk0_26) <- FR_wk0_26$sample_name

# Extracting the desired columns from counts table
selected_counts <- data %>%
  select(all_of(FR_wk0_26$sample_name))

# Print the result
FR_wk0_26$sample_name <- NULL

data_matrix <- as.matrix(round(selected_counts, digits=0))
# Import data into DESEq2 object and define your design
ddsMat <- DESeqDataSetFromMatrix(countData = data_matrix, colData = FR_wk0_26, 
                                 design=~severity)

keep <- rowSums(counts(ddsMat) >= 10) >= 3
ddsMat <- ddsMat[keep,]

geoMean <- apply(counts(ddsMat), 1, gm_mean)
ddsMat <- estimateSizeFactors(ddsMat, geoMean = geoMean)

rld <- rlog(ddsMat, blind=FALSE, fitType = 'local')

makeHCheatmap(rld, FR_wk0_26, severity)
makePCoA(rld, FR_wk0_26, severity)

ddsMat <- DESeq(ddsMat, fitType='local')
resultsNames(ddsMat)
res <- results(ddsMat, contrast = c("severity", "remission", "flare"))
sign_res <- subset(res, padj < 0.05)
nrow(sign_res)

# Subset 3
##################################
# Extracting Patients that were in remission at wk26 
# and turned to flare at wk52 (RF_wk26_52)
joined_data <- new_metadata %>%
  filter(timePoint == "wk26" & severity == "remission") %>%
  inner_join(new_metadata %>% filter(timePoint == "wk52" & severity == "flare"), by = "sample_name")

# Splitting data for wk0 and wk26 and formatting the sample names
wk26_data <- joined_data %>%
  transmute(sample_name = paste0("W26_P", sample_name),
            severity = severity.x)

wk52_data <- joined_data %>%
  transmute(sample_name = paste0("X", sample_name, "_WK52"),
            severity = severity.y)

# Binding rows
RF_wk26_52 <- bind_rows(wk26_data, wk52_data)
row.names(RF_wk26_52) <- RF_wk26_52$sample_name

# Extracting the desired columns from counts table
selected_counts <- data %>%
  select(all_of(RF_wk26_52$sample_name))

# Print the result
RF_wk26_52$sample_name <- NULL

data_matrix <- as.matrix(round(selected_counts, digits=0))
# Import data into DESEq2 object and define your design
ddsMat <- DESeqDataSetFromMatrix(countData = data_matrix, colData = RF_wk26_52, 
                                 design=~severity)

keep <- rowSums(counts(ddsMat) >= 10) >= 3
ddsMat <- ddsMat[keep,]

geoMean <- apply(counts(ddsMat), 1, gm_mean)
ddsMat <- estimateSizeFactors(ddsMat, geoMean = geoMean)

rld <- rlog(ddsMat, blind=FALSE, fitType = 'local')

makeHCheatmap(rld, RF_wk26_52, severity)
makePCoA(rld, RF_wk26_52, severity)

ddsMat <- DESeq(ddsMat, fitType='local')
resultsNames(ddsMat)
res <- results(ddsMat, contrast = c("severity", "flare", "remission"))
sign_res <- subset(res, padj < 0.05)
nrow(sign_res)


# Subset 4
##################################
# Extracting Patients that were in flare at wk26 
# and turned to remission at wk52 (FR_wk26_52)
joined_data <- new_metadata %>%
  filter(timePoint == "wk26" & severity == "flare") %>%
  inner_join(new_metadata %>% filter(timePoint == "wk52" & severity == "remission"), by = "sample_name")

# Splitting data for wk0 and wk26 and formatting the sample names
wk26_data <- joined_data %>%
  transmute(sample_name = paste0("W26_P", sample_name),
            severity = severity.x)

wk52_data <- joined_data %>%
  transmute(sample_name = paste0("X", sample_name, "_WK52"),
            severity = severity.y)

# Binding rows
FR_wk26_52 <- bind_rows(wk26_data, wk52_data)
row.names(FR_wk26_52) <- FR_wk26_52$sample_name

# Extracting the desired columns from counts table
selected_counts <- data %>%
  select(all_of(FR_wk26_52$sample_name))

# Print the result
FR_wk26_52$sample_name <- NULL

data_matrix <- as.matrix(round(selected_counts, digits=0))
# Import data into DESEq2 object and define your design
ddsMat <- DESeqDataSetFromMatrix(countData = data_matrix, 
                                 colData = FR_wk26_52, 
                                 design=~severity)

keep <- rowSums(counts(ddsMat) >= 10) >= 3
ddsMat <- ddsMat[keep,]

geoMean <- apply(counts(ddsMat), 1, gm_mean)
ddsMat <- estimateSizeFactors(ddsMat, geoMean = geoMean)

rld <- rlog(ddsMat, blind=FALSE, fitType = 'local')

makeHCheatmap(rld, FR_wk26_52, severity)
makePCoA(rld, FR_wk26_52, severity)

ddsMat <- DESeq(ddsMat, fitType='local')
resultsNames(ddsMat)
res <- results(ddsMat, contrast = c("severity", "remission","flare"))
sign_res <- subset(res, padj < 0.05)
nrow(sign_res)


# Subset 5
##################################
# Extracting Patients that were in remission at wk0 
# and turned to flare at wk52 (RF_wk0_52)
joined_data <- new_metadata %>%
  filter(timePoint == "wk0" & severity == "remission") %>%
  inner_join(new_metadata %>% filter(timePoint == "wk52" & severity == "flare"), by = "sample_name")

# Splitting data for wk0 and wk26 and formatting the sample names
wk0_data <- joined_data %>%
  transmute(sample_name = paste0("X", sample_name, "_Baseline"),
            severity = severity.x)

wk52_data <- joined_data %>%
  transmute(sample_name = paste0("X", sample_name, "_WK52"),
            severity = severity.y)

# Binding rows
RF_wk0_52 <- bind_rows(wk0_data, wk52_data)
row.names(RF_wk0_52) <- RF_wk0_52$sample_name

# Extracting the desired columns from counts table
selected_counts <- data %>%
  select(all_of(RF_wk0_52$sample_name))

# Print the result
RF_wk0_52$sample_name <- NULL

data_matrix <- as.matrix(round(selected_counts, digits=0))
# Import data into DESEq2 object and define your design
ddsMat <- DESeqDataSetFromMatrix(countData = data_matrix, colData = RF_wk0_52, 
                                 design=~severity)

keep <- rowSums(counts(ddsMat) >= 10) >= 3
ddsMat <- ddsMat[keep,]

# Transform data using the harmonic mean
geoMean <- apply(counts(ddsMat), 1, gm_mean)
ddsMat  <- estimateSizeFactors(ddsMat, geoMean = geoMean)

rld <- rlog(ddsMat, blind=FALSE, fitType = 'local')

makeHCheatmap(rld, RF_wk0_52, severity)
makePCoA(rld, RF_wk0_52, severity)

ddsMat <- DESeq(ddsMat, fitType='local')
resultsNames(ddsMat)
res <- results(ddsMat, contrast = c("severity", "flare", "remission"))
sign_res <- subset(res, padj < 0.05)
nrow(sign_res)


# Subset 6
##################################
# Extracting Patients that were in flare at wk0 
# and turned to remission at wk52 (FR_wk0_52)
joined_data <- new_metadata %>%
  filter(timePoint == "wk0" & severity == "remission") %>%
  inner_join(new_metadata %>% filter(timePoint == "wk52" & severity == "flare"), by = "sample_name")

# Splitting data for wk0 and wk26 and formatting the sample names
wk0_data <- joined_data %>%
  transmute(sample_name = paste0("X", sample_name, "_Baseline"),
            severity = severity.x)

wk52_data <- joined_data %>%
  transmute(sample_name = paste0("X", sample_name, "_WK52"),
            severity = severity.y)

# Binding rows
FR_wk0_52 <- bind_rows(wk0_data, wk52_data)
row.names(FR_wk0_52) <- FR_wk0_52$sample_name

# Extracting the desired columns from counts table
selected_counts <- data %>%
  select(all_of(FR_wk0_52$sample_name))

# Print the result
FR_wk0_52$sample_name <- NULL

data_matrix <- as.matrix(round(selected_counts, digits=0))
# Import data into DESEq2 object and define your design
ddsMat <- DESeqDataSetFromMatrix(countData = data_matrix, colData = RF_wk0_52, 
                                 design=~severity)

keep <- rowSums(counts(ddsMat) >= 10) >= 3
ddsMat <- ddsMat[keep,]

# Transform data using the harmonic mean
geoMean <- apply(counts(ddsMat), 1, gm_mean)
ddsMat  <- estimateSizeFactors(ddsMat, geoMean = geoMean)

rld <- rlog(ddsMat, blind=FALSE, fitType = 'local')

makeHCheatmap(rld, RF_wk0_52, severity)
makePCoA(rld, RF_wk0_52, severity)

ddsMat <- DESeq(ddsMat, fitType='local')
resultsNames(ddsMat)
res <- results(ddsMat, contrast = c("severity", "remission", "flare"))
sign_res <- subset(res, padj < 0.05)
nrow(sign_res)
