# ============================================
# Load R packages 
# ============================================
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)

# ============================================
# File preparation
# ============================================

promoter_window <- fread("/results/04_promoter_accessible_genes/heacon5_promoter_window.bed")
colnames(promoter_window) <- c("chr", "start", "end", "transcript", "offset", "direction")
summit <- fread("/scripts/04_promoter_accessible_genes/heacon5_IDR_filtered_peaks_summit.tsv", sep = "\t")

# ============================================
# ATAC-seq peak summits across all promoter windows and plot relative to the transcription start site
# ============================================

# Define function to calculate peaks summit counts relative to window for each TSS
get_relative_positions <- function(row_as_list, all_summits_df) {
  # Filter summits by chromosome
  summits_on_chr_df <- all_summits_df[all_summits_df$chr == row_as_list$chr, ]
  summits <- summits_on_chr_df$summit
  
  # Define TSS position (relative zero point)
  # '+' strand: TSS (relative 0) is at end - 500 (changed from 100)
  # '-' strand: TSS (relative 0) is at start + 500 (changed from 100)
  
  # Core window: summits between row_as_list$start and row_as_list$end.
  # With end-start = 1500 and TSS offset 500, this window now naturally spans -1000 to +500 relative.
  summits_for_core_processing <- summits[summits >= row_as_list$start & summits <= row_as_list$end]
  
  # Flank window: (start-2000) to (end+2000) - genomic definition remains the same
  flank_start_genomic <- max(0, row_as_list$start - 2000) 
  flank_end_genomic <- row_as_list$end + 2000 # Genomic end of flank
  summits_for_flank_processing <- summits[summits >= flank_start_genomic & summits <= flank_end_genomic]
  
  if (row_as_list$direction == '+') {
    tss_zero_point <- row_as_list$end - 500 # Changed offset to 500
    
    rel_pos_core <- summits_for_core_processing - tss_zero_point
    rel_pos_flank <- summits_for_flank_processing - tss_zero_point
  } else { # direction == '-'
    tss_zero_point <- row_as_list$start + 500 # Changed offset to 500
    
    rel_pos_core <- tss_zero_point - summits_for_core_processing
    rel_pos_flank <- tss_zero_point - summits_for_flank_processing
  }
  
  return(list(core = rel_pos_core, flank = rel_pos_flank))
}

# Apply the function to each row of promoter_window
rel_data_accumulator <- vector("list", nrow(promoter_window)) 

for (i in 1:nrow(promoter_window)) {
  current_row_as_list <- as.list(promoter_window[i, ])
  rel_data_accumulator[[i]] <- get_relative_positions(current_row_as_list, summit)
}

# Flatten and count frequencies for visualization
all_core_positions <- unlist(lapply(rel_data_accumulator, function(x) x$core))
all_flank_positions <- unlist(lapply(rel_data_accumulator, function(x) x$flank))

# Define x-axis range for histogramming (bin centers)
# Range remains -3000 to +2500 as per previous script
x_bins_values <- seq(-3000, 2500, by = 1) 

# Calculate frequency counts for these bins
core_counts_table <- table(factor(all_core_positions, levels = x_bins_values))
flank_counts_table <- table(factor(all_flank_positions, levels = x_bins_values))

# Convert tables to numeric vectors for plotting (these are the y-values)
core_hist_y_values <- as.numeric(core_counts_table)
flank_hist_y_values <- as.numeric(flank_counts_table)

# The x-values for the plot are the bin centers themselves
plot_x_values <- x_bins_values 

# Plot
plot(plot_x_values,flank_hist_y_values, 
     type = 'l', 
     col = rgb(0.5, 0.5, 0.5, alpha = 0.6), 
     xlab = 'Relative Position to TSS (0 = Defined TSS)', 
     ylab = 'Peak Count',
     main = 'Peak Count Across Promoter',
     ylim = range(c(0, core_hist_y_values, flank_hist_y_values)), 
     xlim = c(-3000, 2500), # Plot range remains -3000 to +2500
     xaxt = 'n' 
)
axis(1, at = pretty(plot_x_values, n=10), labels = pretty(plot_x_values, n=10)) 

lines(plot_x_values, core_hist_y_values, col = 'red')
abline(v = 0, col = 'black', lty = 2) 

loess_fit <- loess(flank_hist_y_values ~ plot_x_values,span = 0.1)  # Fit loess
loess_y <- predict(loess_fit)                            # Predicted smooth values
lines(plot_x_values, loess_y, col = 'blue', lwd = 2)     # Add smooth line to plot



# Legend text remains consistent with the relative ranges achieved.
# Core TSS Window is now naturally -1kb to +0.5kb relative.
# Flank plotted range is -3kb to +2.5kb relative.
legend("topright", 
       legend = c('Flank Region', 'Promoter (-1kb to +0.5kb)'), 
       col = c(rgb(0.5, 0.5, 0.5, alpha = 0.6), 'red'), 
       lty = 1, 
       cex = 0.8 
)
grid()


# ============================================
# plot the proportion of peaks within promoter region across different peak sets.
# ============================================

# Read peak distribution summary from HOMER annotation results
# (generated by nf-core/atacseq pipeline)
peak_count_df <- read.table(
  text = "
SampleID   TSS    TTS   Exon  Introns  Intergenic
HCA1      18658  10474 15896  47117     41699
HCA2      25161  14662 26268  70202     62791
HCA       24609  14357 24732  67159     60155
IDR        4002   1631  2438   7287      7564
",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Calculate total and non-TSS peak counts for each sample
peak_count_df$total <- rowSums(peak_count_df[, 2:6])
peak_count_df$nonTSS <- peak_count_df$total - peak_count_df$TSS

# Perform Fisher’s exact test to compare TSS enrichment
# between IDR peaks and the mean of non-IDR samples
fisher_matrix <- matrix(
  c(
    peak_count_df[peak_count_df$SampleID == "IDR", "TSS"],
    peak_count_df[peak_count_df$SampleID == "IDR", "nonTSS"],
    mean(peak_count_df[peak_count_df$SampleID != "IDR", "TSS"]),
    mean(peak_count_df[peak_count_df$SampleID != "IDR", "nonTSS"])
  ),
  nrow = 2,
  byrow = TRUE
)
colnames(fisher_matrix) <- c("TSS", "nonTSS")
rownames(fisher_matrix) <- c("IDR", "nonIDR")

# Run Fisher’s exact test
fisher.test(fisher_matrix)

# Prepare data for plotting
plot_df <- peak_count_df %>%
  rowwise() %>%
  mutate(
    total = sum(c_across(TSS:Intergenic)),
    proportion = TSS / total,
    signif = ifelse(SampleID == "IDR", "*", "")
  )

plot_df$SampleID <- factor(plot_df$SampleID, levels = c("HCA1", "HCA2", "HCA", "IDR"))

# Plot TSS peak proportions across samples
ggplot(plot_df, aes(x = SampleID, y = proportion, group = 1)) +
  geom_line(color = "black") +
  geom_point(size = 3, color = "black") +
  # geom_text(aes(label = signif),
  #           vjust = -2.5, size = 6, fontface = "bold") +
  ylab("Proportion of peaks in TSS") +
  xlab("Sample") +
  coord_cartesian(ylim = c(0.12, 0.18)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# ============================================
# Distribution of biological replicates peaks (HCA1, HCA2), merged peaks (HCA) and reproducible peaks (IDR) across genomic features.
# ============================================

# Define the desired feature order (consistent with final plot)
categories <- c("Intergenic", "TTS", "intron", "exon", "Promoter")

# Clean and simplify annotation column
clean_ann <- function(x) {
  # Define Promoter region as -1000 to +500 from TSS
  x$Feature_type <- ifelse(x$Distance.to.TSS >= -1001 & x$Distance.to.TSS <= 501,
                           "Promoter", x$Annotation)
  # Remove text inside parentheses, e.g. "Exon (1 of many)" → "Exon"
  x$Feature_type <- gsub(" \\(.*\\)", "", x$Feature_type)
  # Remove dashes and suffixes, e.g. "intron -2kb" → "intron"
  x$Feature_type <- gsub(" ?-.*", "", x$Feature_type)
  # Standardize capitalization
  x$Feature_type[x$Feature_type == "promoter"] <- "Promoter"
  x
}

# Read and clean annotation files for four samples
HCA  <- clean_ann(read.delim("/scripts/04_promoter_accessible_genes/heacon5_HCA_annotatePeaks.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE))
HCA1 <- clean_ann(read.delim("/scripts/04_promoter_accessible_genes/heacon5_HCA1_annotatePeaks.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE))
HCA2 <- clean_ann(read.delim("/scripts/04_promoter_accessible_genes/heacon5_HCA2_annotatePeaks.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE))
IDR  <- clean_ann(read.delim("/scripts/04_promoter_accessible_genes/heacon5_IDR_filtered_annotatePeaks.txt",  header = TRUE, sep = "\t", stringsAsFactors = FALSE))

# Count features per category and fill missing ones with 0
count_categories <- function(tbl, cats) {
  cnt <- table(factor(tbl$Feature_type, levels = cats))
  as.integer(cnt)
}

# Combine counts for all samples into a wide table
df_annotation_t <- data.frame(
  Sample      = c("HCA", "HCA1", "HCA2", "IDR"),
  Intergenic  = c(count_categories(HCA,  categories)[1],
                  count_categories(HCA1, categories)[1],
                  count_categories(HCA2, categories)[1],
                  count_categories(IDR,  categories)[1]),
  TTS         = c(count_categories(HCA,  categories)[2],
                  count_categories(HCA1, categories)[2],
                  count_categories(HCA2, categories)[2],
                  count_categories(IDR,  categories)[2]),
  intron      = c(count_categories(HCA,  categories)[3],
                  count_categories(HCA1, categories)[3],
                  count_categories(HCA2, categories)[3],
                  count_categories(IDR,  categories)[3]),
  exon        = c(count_categories(HCA,  categories)[4],
                  count_categories(HCA1, categories)[4],
                  count_categories(HCA2, categories)[4],
                  count_categories(IDR,  categories)[4]),
  Promoter    = c(count_categories(HCA,  categories)[5],
                  count_categories(HCA1, categories)[5],
                  count_categories(HCA2, categories)[5],
                  count_categories(IDR,  categories)[5]),
  check.names = FALSE
)

# Transform to long format for ggplot visualisation
df_long <- df_annotation_t %>%
  pivot_longer(cols = -Sample, names_to = "Feature", values_to = "Count") %>%
  mutate(
    Sample  = factor(Sample,  levels = c("IDR","HCA","HCA2","HCA1")),
    Feature = factor(Feature, levels = categories)
  )

# Plot stacked bar chart (percentage + count labels, horizontal layout)
ggplot(df_long, aes(x = Sample, y = Count, fill = Feature)) +
  geom_bar(stat = "identity", position = "fill") +  # normalize to proportions
  geom_text(aes(label = Count),
            position = position_fill(vjust = 0.5),
            size = 3, color = "white") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Proportion",
    title = "Peak Annotation Proportion Across Samples"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_blank()
  )


