library(effectsize)  # Ensure this library is loaded for Cohen's d

setwd('/Users/juanjovel/jj/data_analysis/heatherArmstrong/living_w_IBD/mixOmics/analyses_230913')

datafile <- 'taxa_only_0-52w_100up_CPM.tsv'
metadatafile <- 'meta.tsv'

data <- read.table(datafile, sep = '\t', row.names = 1, header = T)
metadata <- read.table(metadatafile, sep = '\t', row.names = 1, header = T)

# Initialize new_df with the appropriate structure
new_df <- data.frame(matrix(ncol = length(row.names(metadata)), nrow = nrow(data)))
colnames(new_df) <- row.names(metadata)  # Name columns after samples
row.names(new_df) <- row.names(data)

# Loop over each sample
for (sample in row.names(metadata)) {
  sample0_colname <- paste0(sample, "0")  # Create name for baseline column
  sample52_colname <- paste0(sample, "52")  # Create name for week 52 column
  
  if (sample0_colname %in% colnames(data) && sample52_colname %in% colnames(data)) {
    new_df[, sample] <- data[, sample52_colname] - data[, sample0_colname]
  } else {
    warning(paste("Columns", sample0_colname, "or", sample52_colname, "not found in data"))
  }
}


# Convert metadata to a vector with names for easy grouping

# Assuming the metadata and new_df are correctly loaded and formatted
group_info <- metadata$group
names(group_info) <- row.names(metadata)

# Convert all data in new_df to numeric
new_df[] <- lapply(new_df, function(x) as.numeric(as.character(x)))

write.table(new_df, "diff_table.txt", sep = '\t', quote = F)

# Prepare a dataframe to store results
# Adding a simple filtering criterion: taxa must be present at a minimum level in at least a few samples
# Parameters for filtering
min_samples_presence <- 3      # Minimum number of samples that must exceed the threshold
min_presence_threshold <- 0.1  # The minimum value to be considered as present

# Function to determine if a row should be kept
should_keep_row <- function(row) {
  sum(row > min_presence_threshold) >= min_samples_presence
}

# Apply this function to each row of new_df to create a logical index of rows to keep
rows_to_keep <- apply(new_df, 1, should_keep_row)

# Filter the dataframe based on this logical index of rows
filtered_new_df <- new_df[rows_to_keep, ]

# Initialize the results dataframe properly
results_df <- data.frame(
  taxa = row.names(filtered_new_df),  # This is correct and should match the number of rows
  p_value = numeric(nrow(filtered_new_df)),  # Number of elements should match the number of rows (taxa)
  effect_size = numeric(nrow(filtered_new_df)),  # Same here, as each row is a taxa, each needs a result slot
  stringsAsFactors = FALSE
)

# Modify the loop to handle potential empty returns from cohens_d
for (i in seq_along(filtered_new_df)) {
  values <- as.numeric(filtered_new_df[, i])
  print(values)
  values_RF <- values[metadata$group == "RF"]
  values_RR <- values[metadata$group == "RR"]
  
  if (length(values_RF) > 1 && length(values_RR) > 1) {
    test_result <- wilcox.test(values_RF, values_RR, exact = FALSE)
    results_df$p_value[i] <- test_result$p.value
    
    # Calculate Cohen's d, ensure it is handled correctly
    if (!is.na(test_result$p.value)) {
      cohen_d_result <- effectsize::cohens_d(values_RF, values_RR, pooled = TRUE)
      results_df$effect_size[i] <- ifelse(length(cohen_d_result$estimate) > 0, cohen_d_result$estimate, NA)
    } else {
      results_df$effect_size[i] <- NA
    }
  } else {
    results_df$p_value[i] <- NA
    results_df$effect_size[i] <- NA
  }
}


# View the results
print(results_df)

write.table(results, "results_Utest.txt", sep = '\t', quote = F, row.names = F)

write.table(new_df, "diff_data.txt", sep = '\t', quote = F, row.names = T)
                    
