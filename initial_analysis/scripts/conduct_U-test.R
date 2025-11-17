# Load packages
library(readr)
library(dplyr)

setwd("/Users/juanjovel/Library/CloudStorage/OneDrive-UniversityofCalgary/jj/UofC/data_analysis/diego_nobrega/sageGrouse")

# Read the data
df <- read_tsv("assembly_stats.tsv")

# Columns to test (edit this according to your header; remove non-numeric columns!)
cols_to_test <- c(
  "Total length [sum]",
  "not_aligned length",
  "# contigs [sum]",
  "not_aligned # contigs",
  "Largest contig [max]",
  "Total aligned length [sum]",
  "# misassemblies [sum]",
  "Misassembled contigs length [sum]",
  "# mismatches per 100 kbp [wmean]",
  "# indels per 100 kbp [wmean]",
  "# N's per 100 kbp [wmean]",
  "Genome fraction (%) [wmean]",
  "Duplication ratio [wmean]",
  "Total length (>= 1000 bp) [sum]",
  "Total length (>= 10000 bp) [sum]",
  "Total length (>= 50000 bp) [sum]",
  "Largest alignment [max]",
  "NGA50 [max]",
  "LGA50 [max]"
)

# Wilcoxon test for each column
test_results <- lapply(cols_to_test, function(column) {
  # Make sure column exists and is numeric
  if (column %in% names(df) && is.numeric(df[[column]])) {
    group1 <- df %>% filter(group == "captivity_E") %>% pull(column)
    group2 <- df %>% filter(group == "wild_I") %>% pull(column)
    wt <- wilcox.test(group1, group2)
    data.frame(
      column = column,
      p_value = wt$p.value
    )
  } else {
    data.frame(
      column = column,
      p_value = NA
    )
  }
})

# Combine results
results_df <- bind_rows(test_results)

# Print results
print(results_df)

# Save results
write_tsv(results_df, "mann_whitney_results.tsv")

library(ggplot2)

# Metrics to plot
sig_metrics <- c(
  "Total length [sum]",
  "not_aligned length",
  "# contigs [sum]",
  "not_aligned # contigs",
  "Largest contig [max]"
)

# Loop to save plots
for (metric in sig_metrics) {
  plot_df <- df %>%
    group_by(group) %>%
    summarize(
      mean_value = mean(.data[[metric]], na.rm = TRUE),
      sd_value = sd(.data[[metric]], na.rm = TRUE)
    )
  
  p <- ggplot(plot_df, aes(x = group, y = mean_value, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.6) +
    geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
                  width = 0.2, position = position_dodge(0.6)) +
    labs(title = paste("Group comparison:", metric),
         y = metric,
         x = "Group") +
    theme_bw() +
    theme(legend.position = "none")

  ggsave(filename = paste0("barplot_", gsub("[^A-Za-z0-9]", "", metric), ".png"), plot = p)
}

