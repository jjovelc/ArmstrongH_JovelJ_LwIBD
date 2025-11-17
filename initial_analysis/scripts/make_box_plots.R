#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(tools)
})

## -------------------- paths --------------------
# setwd('/Users/juanjovel/Library/CloudStorage/OneDrive-UniversityofCalgary/jj/UofC/data_analysis/diego_nobrega/sageGrouse')
counts_file <- "Greater_Sage-grouse_kraken2_250619_counts_sort.tsv"
meta_file   <- "metadata.tsv"
diffab_file <- "differential_abundance_02_capt_vs_03_sil_q0.05_FC2.tsv"

## -------------------- helpers --------------------
parse_last_taxon <- function(taxa_string) {
  if (is.na(taxa_string) || taxa_string == "") return("Unknown")
  parts <- strsplit(taxa_string, "\\|")[[1]]
  trimws(parts[length(parts)])
}

strip_rank_prefixes <- function(x) {
  # remove rank tags commonly seen in taxonomies
  gsub("(^|_)(k__|p__|c__|o__|f__|g__|s__)", "", x, ignore.case = TRUE)
}

normalize_taxon <- function(x) {
  x <- trimws(x)
  x <- strip_rank_prefixes(x)
  x <- gsub("\\s+", " ", x)
  ifelse(nchar(x) == 0, "Unknown", x)
}

safe_fname <- function(s, n = 150) {
  s <- gsub("[^A-Za-z0-9_.-]+", "_", s)
  substr(s, 1, n)
}

## Colors (from your figure)
group_order <- c("03_sil", "02_capt")
group_cols  <- c("03_sil" = "#ff6347", "02_capt" = "#1e90ff")



## -------------------- read inputs --------------------
message("Reading input files...")
counts <- fread(counts_file)
meta   <- fread(meta_file)
# diffab has no header; name its columns
diffab <- fread(diffab_file, header = FALSE)
setnames(diffab, c("taxon_lineage", "logFC", "pvalue", "qvalue"))

## -------------------- extract taxa (first column) --------------------
taxa_list_raw <- unique(diffab$taxon_lineage)
fwrite(data.table(taxon = taxa_list_raw), "extracted_taxa_from_diffab.tsv", sep = "\t")

## -------------------- build parsed+clean names --------------------
# DIFFAB side
diffab[, taxon_last_with_prefix := vapply(taxon_lineage, parse_last_taxon, character(1))]
diffab[, taxon_parsed_clean     := normalize_taxon(taxon_last_with_prefix)]
# Order taxa by q-value (ascending) to assign Feature_1, Feature_2, ...
suppressWarnings(diffab[, qvalue := as.numeric(qvalue)])
diffab[is.na(qvalue), qvalue := Inf]
setorder(diffab, qvalue, taxon_parsed_clean)

## -------------------- RPM normalize counts --------------------
# RPM = (counts / library_size) * 1e6
tax_col     <- names(counts)[1]               # "Taxa"
sample_cols <- setdiff(names(counts), tax_col)

# enforce numeric
counts[, (sample_cols) := lapply(.SD, as.numeric), .SDcols = sample_cols]
lib_sizes <- counts[, lapply(.SD, sum, na.rm = TRUE), .SDcols = sample_cols]
denoms <- as.numeric(lib_sizes[1, ])
denoms[denoms == 0] <- NA

rpm <- copy(counts)
for (i in seq_along(sample_cols)) {
  sc <- sample_cols[i]
  rpm[[sc]] <- (rpm[[sc]] / denoms[i]) * 1e6
}
rpm[is.na(rpm)] <- 0

# counts side: keep BOTH the last token WITH prefix (for titles) and the cleaned version (for matching)
rpm[, taxon_last_with_prefix := vapply(get(tax_col), parse_last_taxon, character(1))]
rpm[, taxon_parsed_clean     := normalize_taxon(taxon_last_with_prefix)]

## -------------------- long & aggregate by parsed_clean --------------------
long_rpm <- melt(
  rpm,
  id.vars = c(tax_col, "taxon_last_with_prefix", "taxon_parsed_clean"),
  measure.vars = sample_cols,
  variable.name = "sample",
  value.name   = "rpm"
)

# collapse any duplicated rows that map to same parsed_clean x sample
long_rpm_agg <- long_rpm[, .(rpm = sum(rpm, na.rm = TRUE)),
                         by = .(taxon_parsed_clean, sample)]

## -------------------- metadata prep --------------------
# Expect columns: sample_id, group  (label optional)
if (!all(c("sample_id","group") %in% names(meta))) {
  stop("metadata.tsv must contain columns: sample_id and group")
}
meta[, group_norm := gsub("\\s+", "_", trimws(group))]
meta[, sample     := as.character(sample_id)]

meta_sub <- meta[group_norm %in% group_order, .(sample, group = factor(group_norm, levels = group_order))]

## -------------------- join RPM + metadata --------------------
dat <- merge(long_rpm_agg, meta_sub, by = "sample", all = FALSE)

# make a display-name map that PRESERVES the rank prefix for titles
diffab_map <- unique(diffab[, .(taxon_parsed_clean, taxon_last_with_prefix)])
counts_map <- unique(rpm[,    .(taxon_parsed_clean, taxon_last_with_prefix)])
display_map <- merge(diffab_map, counts_map, by = "taxon_parsed_clean", all = TRUE,
                     suffixes = c("_diffab", "_counts"))
display_map[, display_with_prefix :=
              ifelse(!is.na(taxon_last_with_prefix_diffab),
                     taxon_last_with_prefix_diffab,
                     taxon_last_with_prefix_counts)]
display_map[is.na(display_with_prefix), display_with_prefix := taxon_parsed_clean]  # fallback
display_map <- display_map[, .(taxon_parsed_clean, display_with_prefix)]

## -------------------- filter to taxa present in diff list --------------------
taxa_target <- unique(diffab$taxon_parsed_clean)
dat_sub <- dat[taxon_parsed_clean %in% taxa_target]

matched_n <- uniqueN(dat_sub$taxon_parsed_clean)
message("Matched (parsed-clean) taxa: ", matched_n, " of ", uniqueN(taxa_target))

## -------------------- output dir --------------------
outdir <- "boxplots"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## -------------------- plotting (RPM) --------------------
# order taxa for plotting by q-value (ascending)
taxa_order <- diffab$taxon_parsed_clean

# iterate and save one PNG per taxon
n_done <- 0L
for (i in seq_along(taxa_order)) {
  tx <- taxa_order[i]
  dt <- dat_sub[taxon_parsed_clean == tx]
  
  if (nrow(dt) == 0) next
  
  # retrieve title with rank prefix (e.g. s__, g__, f__)
  disp <- display_map[match(tx, display_map$taxon_parsed_clean), display_with_prefix]
  if (is.na(disp) || disp == "") disp <- tx
  
  p <- ggplot(dt, aes(x = group, y = rpm, fill = group)) +
    geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.6) +
    geom_jitter(width = 0.12, size = 1.8, alpha = 0.75, shape = 21, color = "black") +
    scale_fill_manual(values = group_cols, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    labs(
      x = "Group",
      y = "Reads per million (RPM)",
      title = sprintf("Feature_%d: %s", i, disp)
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  fname <- file.path(outdir, sprintf("%04d__%s.png", i, safe_fname(disp)))
  ggsave(filename = fname, plot = p, width = 6, height = 4.5, dpi = 200, bg = "white")
  n_done <- n_done + 1L
}

message("Saved ", n_done, " boxplots into ", outdir, "/")

## ---------- Zero-mean check (RPM) & export ----------
# Mean RPM per taxon (parsed_clean) and group
means_long <- dat_sub[, .(avg = mean(rpm, na.rm = TRUE)),
                      by = .(taxon_parsed_clean, group)]

# Wide table with one column per group
means_wide <- data.table::dcast(
  means_long,
  taxon_parsed_clean ~ group,
  value.var = "avg",
  fill = 0
)

# Ensure both groups exist as columns
for (g in group_order) {
  if (!g %in% names(means_wide)) means_wide[[g]] <- 0
}

# Rename to explicit columns and keep order (02_capt, 03_sil)
data.table::setnames(
  means_wide,
  old = c("02_capt", "03_sil"),
  new = c("avg_02_capt", "avg_03_sil"),
  skip_absent = TRUE
)

# Keep taxa where either group's average RPM is exactly zero
zeros <- means_wide[avg_02_capt == 0 | avg_03_sil == 0]

# Attach display name WITH taxonomic rank prefix for the title/ID
zeros <- merge(zeros, display_map, by = "taxon_parsed_clean", all.x = TRUE)
zeros[is.na(display_with_prefix) | display_with_prefix == "", display_with_prefix := taxon_parsed_clean]
data.table::setcolorder(zeros, c("display_with_prefix", "avg_02_capt", "avg_03_sil"))
data.table::setnames(zeros, "display_with_prefix", "taxon")

# Write out the 3-column file (taxon + two averages in RPM)
fwrite(zeros[, .(taxon, avg_02_capt, avg_03_sil)],
       file = "taxa_with_zero_group_mean.tsv", sep = "\t")

message("Wrote taxa_with_zero_group_mean.tsv (RPM).")
