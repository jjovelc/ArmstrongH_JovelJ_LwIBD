# This script takes as input a series of kraken2 (mpa style) generated
# by running: python parse_taxa.py kraken2_mpa-style_report.tsv
# e.g. python parse_taxa.py inflammation_all_samples_kraken2-counts.tsv
# The data will be normalized, and alpha and beta diversity analyses
# conducted, and plots will be generated

# Make sure you have the required libraries, otherwise:
# BiocManager::install("phyloseq")
# install.packages("tidyverse")
# install.packages("remotes")
# remotes::install_github("vegandevs/vegan")
rm(list = ls())
library(phyloseq)
library(tidyverse)
library(vegan)

# Change directory to the directory where you have the input data
setwd("/Users/juanjovel/Library/CloudStorage/OneDrive-UniversityofCalgary/jj/UofC/data_analysis/diego_nobrega/sageGrouse")

# Define the levels you want to iterate over
levels <- 1:7

# Define the taxonomic ranks
taxonomic_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Function to extract terminal leaf from taxonomic string
extract_terminal_leaf <- function(taxa_name) {
  # Split by | and get the last non-empty element
  parts <- strsplit(taxa_name, "\\|")[[1]]
  parts <- parts[parts != "" & !is.na(parts)]
  if (length(parts) > 0) {
    terminal <- parts[length(parts)]
    # Remove prefixes like k__, p__, c__, etc.
    terminal <- gsub("^[kpcofgs]__", "", terminal)
    return(terminal)
  } else {
    return("Unknown")
  }
}

# Function to perform pairwise alpha diversity analysis
perform_pairwise_alpha_analysis <- function(physeq_raw, level_suffix, rank_name_for_file, taxonomic_rank_name) {
  # Define pairwise comparisons
  pairwise_comparisons <- list(
    c("00_wild", "02_capt_ant"),
    c("00_wild", "01_capt_no_ant"),
    c("01_capt_no_ant", "02_capt_ant")
  )
  
  comparison_names <- c("wild_vs_capt_ant", "wild_vs_capt_no_ant", "capt_no_ant_vs_capt_ant")
  
  indices <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
  
  for (i in seq_along(pairwise_comparisons)) {
    group1 <- pairwise_comparisons[[i]][1]
    group2 <- pairwise_comparisons[[i]][2]
    comparison_name <- comparison_names[i]
    
    # Subset phyloseq object to only include the two groups being compared
    sample_data_subset <- sample_data(physeq_raw)
    keep_samples <- rownames(sample_data_subset)[sample_data_subset$Condition %in% c(group1, group2)]
    physeq_subset <- prune_samples(keep_samples, physeq_raw)
    
    # Perform alpha diversity analysis for each index
    for (index in indices) {
      # Create plot
      p <- plot_richness(physeq_subset, x = "Condition", measures = index) + 
        geom_boxplot(aes(color = Condition)) +
        geom_jitter(width = 0.2, alpha = 0.7) +
        ggtitle(paste(level_suffix, index, ":", group1, "vs", group2, sep = " ")) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 12, hjust = 0.5)
        )
      
      # Save plot
      plot_filename <- paste("alpha_diversity_", rank_name_for_file, "_", index, "_", comparison_name, ".png", sep = "")
      ggsave(plot_filename, plot = p, width = 8, height = 6)
      
      # Perform statistical test (Wilcoxon rank-sum test)
      alpha_values <- estimate_richness(physeq_subset, measures = index)
      alpha_values$Condition <- sample_data(physeq_subset)$Condition
      
      group1_values <- alpha_values[alpha_values$Condition == group1, index]
      group2_values <- alpha_values[alpha_values$Condition == group2, index]
      
      # Wilcoxon test
      wilcox_result <- wilcox.test(group1_values, group2_values)
      
      # Print results
      cat(sprintf("\n%s %s %s vs %s - %s Index:\n", level_suffix, taxonomic_rank_name, group1, group2, index))
      cat(sprintf("Group %s (n=%d): Mean = %.3f, Median = %.3f\n", 
                  group1, length(group1_values), mean(group1_values), median(group1_values)))
      cat(sprintf("Group %s (n=%d): Mean = %.3f, Median = %.3f\n", 
                  group2, length(group2_values), mean(group2_values), median(group2_values)))
      cat(sprintf("Wilcoxon test: W = %.3f, p-value = %.4f\n", 
                  wilcox_result$statistic, wilcox_result$p.value))
    }
  }
}

# Function to perform pairwise beta diversity analysis
perform_pairwise_beta_analysis <- function(physeq_relative, level_suffix, rank_name_for_file, taxonomic_rank_name) {
  # Define pairwise comparisons
  pairwise_comparisons <- list(
    c("00_wild", "02_capt_ant"),
    c("00_wild", "01_capt_no_ant"),
    c("01_capt_no_ant", "02_capt_ant")
  )
  
  comparison_names <- c("wild_vs_capt_ant", "wild_vs_capt_no_ant", "capt_no_ant_vs_capt_ant")
  
  distances <- c("bray", "jaccard", "jsd")
  
  for (i in seq_along(pairwise_comparisons)) {
    group1 <- pairwise_comparisons[[i]][1]
    group2 <- pairwise_comparisons[[i]][2]
    comparison_name <- comparison_names[i]
    
    # Subset phyloseq object to only include the two groups being compared
    sample_data_subset <- sample_data(physeq_relative)
    keep_samples <- rownames(sample_data_subset)[sample_data_subset$Condition %in% c(group1, group2)]
    physeq_subset <- prune_samples(keep_samples, physeq_relative)
    
    # Get sample metadata for subset
    sample_metadata_subset <- sample_data(physeq_subset)
    
    for (dist in distances) {
      # Calculate distance matrix for subset
      distance_matrix <- phyloseq::distance(physeq_subset, method = dist)
      
      # Perform PCoA
      pcoa_results <- cmdscale(distance_matrix, eig = TRUE)
      pcoa_df <- as.data.frame(pcoa_results$points)
      pcoa_df$condition <- sample_metadata_subset$Condition
      pcoa_df$sample_id <- rownames(pcoa_df)
      
      # Create proper data frame for PERMANOVA
      sample_metadata_df_subset <- data.frame(
        SampleID = rownames(sample_metadata_subset),
        Condition = sample_metadata_subset$Condition,
        row.names = rownames(sample_metadata_subset)
      )
      
      # Debug: print the structure of the data frame
      cat(sprintf("Debug: sample_metadata_df_subset structure for %s vs %s:\n", group1, group2))
      cat(sprintf("Class: %s\n", class(sample_metadata_df_subset)))
      cat(sprintf("Dimensions: %d x %d\n", nrow(sample_metadata_df_subset), ncol(sample_metadata_df_subset)))
      cat(sprintf("Column names: %s\n", paste(colnames(sample_metadata_df_subset), collapse = ", ")))
      cat(sprintf("Condition values: %s\n", paste(sample_metadata_df_subset$Condition, collapse = ", ")))
      
      # Conduct PERMANOVA
      permanova_results <- adonis2(distance_matrix ~ Condition, data = sample_metadata_df_subset)
      
      # Conduct ANOSIM
      anosim_results <- anosim(distance_matrix, sample_metadata_subset$Condition)
      
      # Extract p-value from PERMANOVA results
      permanova_pvalue <- permanova_results$`Pr(>F)`[1]
      
      # Add error checking for p-value extraction
      if (is.null(permanova_pvalue) || is.na(permanova_pvalue)) {
        permanova_pvalue <- "NA"
        cat("Warning: Could not extract PERMANOVA p-value\n")
      }
      
      # Print results
      cat(sprintf("\n%s %s %s vs %s - %s Distance:\n", level_suffix, taxonomic_rank_name, group1, group2, dist))
      cat(sprintf("PERMANOVA Results:\n"))
      print(permanova_results)
      cat(sprintf("ANOSIM Results:\n"))
      print(anosim_results)
      
      # Create PCoA plot
      p <- ggplot(pcoa_df, aes(x = V1, y = V2, color = condition)) +
        geom_point(size = 4) +
        geom_text(aes(label = sample_id), nudge_x = 0.02, nudge_y = 0.02, 
                  check_overlap = TRUE, size = 3) +
        theme_bw() +
        theme(
          legend.position = 'top', 
          legend.title = element_blank(),
          plot.title = element_text(size = 12, hjust = 0.5)
        ) +
        ggtitle(paste(level_suffix, dist, ":", group1, "vs", group2, 
                     "\nPERMANOVA p =", format(permanova_pvalue, digits = 4),
                     "ANOSIM R =", format(anosim_results$statistic, digits = 4)))
      
      # Save plot
      plot_filename <- paste("pcoa_", rank_name_for_file, "_", dist, "_", comparison_name, ".png", sep = "")
      ggsave(plot_filename, plot = p, width = 10, height = 8)
    }
  }
}

# Iterate over each level
for (level in levels) {
  # Generate file name based on level
  
  infile <- sprintf("Greater_Sage-grouse_kraken2_250619_counts_reordered_level_%d.tsv", level)
  
  # Extract level suffix and number
  level_suffix <- sprintf('level_%d', level)
  taxonomic_rank_name <- taxonomic_ranks[level]
  rank_name_for_file <- gsub(" ", "_", tolower(taxonomic_rank_name))
  
  # Read the OTU table from the file
  shotgun_table <- read.csv(infile, sep = '\t', header = TRUE, check.names = FALSE, row.names = 1)
  
  # Ensure taxa names in OTU table are consistent and formatted correctly
  taxa_names_in_otu_table <- rownames(shotgun_table)
  
  # For diversity indices calculations requiring raw counts
  otu_table_obj_raw <- otu_table(as.matrix(shotgun_table), taxa_are_rows = TRUE)
  
  # For parts of analysis where relative abundances are appropriate
  shotgun_table_relative_abundance <- apply(shotgun_table, 2, function(x) x / sum(x))
  otu_table_obj_relative <- otu_table(as.matrix(shotgun_table_relative_abundance), taxa_are_rows = TRUE)
  
  # Select the taxonomic ranks to use based on the level
  ranks_to_use <- taxonomic_ranks[1:level]
  
  # Create taxa information
  taxa_info <- data.frame(Taxa = taxa_names_in_otu_table) %>%
    separate(Taxa, into = ranks_to_use, sep = "\\|", fill = "right") %>%
    as.data.frame()
  
  row.names(taxa_info) <- taxa_names_in_otu_table
  
  tax_table_obj <- tax_table(as.matrix(taxa_info))
  
  samples <- colnames(shotgun_table)
  
  # Metadata
  sample_metadata_df <- data.frame(
    SampleID = samples,
  
    Condition = c("00_wild","00_wild","00_wild","00_wild","00_wild","02_capt_ant","02_capt_ant","02_capt_ant","02_capt_ant",
                  "02_capt_ant","02_capt_ant","02_capt_ant","02_capt_ant","02_capt_ant","02_capt_ant","02_capt_ant","02_capt_ant",
                  "02_capt_ant","02_capt_ant","02_capt_ant","01_capt_no_ant","01_capt_no_ant","01_capt_no_ant","01_capt_no_ant",
                  "01_capt_no_ant","01_capt_no_ant","01_capt_no_ant","01_capt_no_ant","01_capt_no_ant","01_capt_no_ant")
    
  )
  
  row.names(sample_metadata_df) <- sample_metadata_df$SampleID
  
  sample_metadata <- sample_data(sample_metadata_df)
  
  # Create phyloseq object for raw counts
  physeq_raw <- phyloseq(otu_table_obj_raw, tax_table_obj, sample_metadata)
  
  # Create phyloseq object for relative abundances
  # This object can be used for beta diversity analysis or other analyses where relative abundances are appropriate
  physeq_relative <- phyloseq(otu_table_obj_relative, tax_table_obj, sample_metadata)
  
  # Generate stacked bar plots for levels 2-7 (Phylum through Species)
  if (level >= 2) {
    cat(sprintf("\nGenerating stacked bar plot for %s (Level %d)...\n", taxonomic_rank_name, level))
    
    # Calculate mean relative abundance across all samples for each taxon
    mean_abundance <- rowMeans(shotgun_table_relative_abundance)
    
    # Get top 20 most abundant taxa (or max available if less than 20)
    n_taxa <- min(20, length(mean_abundance))
    top_taxa_indices <- order(mean_abundance, decreasing = TRUE)[1:n_taxa]
    top_taxa_names <- names(mean_abundance)[top_taxa_indices]
    
    # Subset data to top taxa
    top_taxa_data <- shotgun_table_relative_abundance[top_taxa_names, , drop = FALSE]
    
    # Extract terminal leaves for legend
    terminal_leaves <- sapply(rownames(top_taxa_data), extract_terminal_leaf)
    
    # Create data frame for plotting
    plot_data <- top_taxa_data %>%
      as.data.frame() %>%
      rownames_to_column("Taxa") %>%
      mutate(Terminal_Leaf = terminal_leaves) %>%
      pivot_longer(cols = -c(Taxa, Terminal_Leaf), names_to = "Sample", values_to = "Relative_Abundance") %>%
      left_join(sample_metadata_df, by = c("Sample" = "SampleID"))
    
    # Create discrete color palette with support for >20 colors
    n_colors <- length(unique(plot_data$Terminal_Leaf))
    
    # Define a comprehensive palette combining multiple discrete palettes
    discrete_palette <- c(
      # RColorBrewer Set3 (12 colors)
      RColorBrewer::brewer.pal(12, "Set3"),
      # RColorBrewer Paired (12 colors)
      RColorBrewer::brewer.pal(12, "Paired"),
      # RColorBrewer Set1 (9 colors)
      RColorBrewer::brewer.pal(9, "Set1"),
      # RColorBrewer Set2 (8 colors)  
      RColorBrewer::brewer.pal(8, "Set2"),
      # RColorBrewer Dark2 (8 colors)
      RColorBrewer::brewer.pal(8, "Dark2"),
      # RColorBrewer Accent (8 colors)
      RColorBrewer::brewer.pal(8, "Accent"),
      # Additional pastel colors
      "#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF", "#E6BAFF",
      "#FFCCCB", "#FFDAB9", "#E0E0E0", "#B0E0E6", "#F0E68C", "#DDA0DD",
      "#98FB98", "#F5DEB3", "#D3D3D3", "#AFEEEE", "#FAFAD2", "#THISTLE"
    )
    
    # Remove duplicates and select needed colors
    discrete_palette <- unique(discrete_palette)
    colors <- discrete_palette[1:min(n_colors, length(discrete_palette))]
    
    # If we still need more colors, supplement with generated pastels
    if (n_colors > length(discrete_palette)) {
      additional_needed <- n_colors - length(discrete_palette)
      hues <- seq(0, 1, length.out = additional_needed + 1)[1:additional_needed]
      additional_colors <- hsv(h = hues, s = 0.35, v = 0.85)
      colors <- c(colors, additional_colors)
    }
    
    # Create stacked bar plot
    p_stacked <- ggplot(plot_data, aes(x = Sample, y = Relative_Abundance, fill = Terminal_Leaf)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = colors) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, hjust = 0.5)
      ) +
      labs(
        x = "Sample",
        y = "Relative Abundance",
        fill = taxonomic_rank_name,
        title = paste("Top", n_taxa, "Most Abundant", taxonomic_rank_name, "- Stacked Bar Plot")
      ) +
      facet_wrap(~ Condition, scales = "free_x", nrow = 1)
    
    # Save stacked bar plot
    stacked_plot_filename <- paste("stacked_barplot_", rank_name_for_file, "_top", n_taxa, ".png", sep = "")
    ggsave(stacked_plot_filename, plot = p_stacked, width = 12, height = 8, dpi = 300)
    
    cat(sprintf("Stacked bar plot saved as: %s\n", stacked_plot_filename))
  }
  
  # Overall alpha diversity calculations (all groups together)
  cat(sprintf("\nPerforming overall alpha diversity analysis for %s (Level %d)...\n", taxonomic_rank_name, level))
  indices <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
  for (index in indices) {
    p <- plot_richness(physeq_raw, x = "Condition", measures = index) + 
      geom_boxplot(aes(color = Condition)) +
      ggtitle(paste(level_suffix, index, sep = "_"))
    
    # Save plot
    plot_filename <- paste("alpha_diversity_", rank_name_for_file, "_", index, ".png", sep = "")
    ggsave(plot_filename, plot = p, width = 10, height = 8)
  }
  
  # Pairwise alpha diversity analysis
  cat(sprintf("\nPerforming pairwise alpha diversity analysis for %s (Level %d)...\n", taxonomic_rank_name, level))
  perform_pairwise_alpha_analysis(physeq_raw, level_suffix, rank_name_for_file, taxonomic_rank_name)
  
  # Overall beta diversity analysis (all groups together)
  cat(sprintf("\nPerforming overall beta diversity analysis for %s (Level %d)...\n", taxonomic_rank_name, level))
  distances <- c("bray", "jaccard", "jsd")
  for (dist in distances) {
    distance_matrix <- phyloseq::distance(physeq_relative, method = dist)
    pcoa_results <- cmdscale(distance_matrix, eig = TRUE)
    pcoa_df <- as.data.frame(pcoa_results$points)
    pcoa_df$condition <- sample_metadata_df$Condition
    
    # Conduct PERMANOVA using adonis2
    permanova_results <- adonis2(distance_matrix ~ Condition, data = sample_metadata_df)
    # Conduct ANOSIM
    anosim_results <- anosim(distance_matrix, sample_metadata_df$Condition)
    
    # Extract p-value from PERMANOVA results
    permanova_pvalue <- permanova_results$`Pr(>F)`[1]
    
    # Add error checking for p-value extraction
    if (is.null(permanova_pvalue) || is.na(permanova_pvalue)) {
      permanova_pvalue <- "NA"
      cat("Warning: Could not extract PERMANOVA p-value\n")
    }
    
    # Print results
    cat(sprintf("\n%s %s PERMANOVA Results:\n", level_suffix, dist))
    print(permanova_results)
    cat(sprintf("\n%s %s ANOSIM Results:\n", level_suffix, dist))
    print(anosim_results)
    
    # Plot PCoA with corrected p-value access
    p <- ggplot(pcoa_df, aes(x = V1, y = V2, color = condition)) +
      geom_point(size = 5) +
      geom_text(aes(label = row.names(pcoa_df)), nudge_x = 0.02, nudge_y = 0.02, check_overlap = TRUE) +
      theme_bw() +
      theme(legend.position = 'top', legend.title = element_blank()) +
      ggtitle(paste(level_suffix, dist, "PERMANOVA p =", format(permanova_pvalue, digits = 4), "ANOSIM R =", format(anosim_results$statistic, digits = 4)))
    
    # Save plot with results
    plot_filename <- paste("pcoa_", rank_name_for_file, "_", dist, "_with_stats.png", sep = "")
    ggsave(plot_filename, plot = p, width = 10, height = 8)
  }
  
  # Pairwise beta diversity analysis
  cat(sprintf("\nPerforming pairwise beta diversity analysis for %s (Level %d)...\n", taxonomic_rank_name, level))
  perform_pairwise_beta_analysis(physeq_relative, level_suffix, rank_name_for_file, taxonomic_rank_name)
}
