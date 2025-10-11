# ============================================
# Load R packages 
# ============================================
library(dplyr)
library(readxl)
library(ggplot2)
library(scales)

# ============================================
# File preparation
# ============================================

# read data frame from last step
ess_ml <- read.table("/results/06_gene_essentiality_analysis/ess_ml.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# read genes modules file 
ise_gene_modules <- read.csv(
  "/scripts/07_accessibility_in_coexpression_modules/ise_genes_modules.csv",
  header = FALSE, col.names = c("ise", "cluster")
)

# add modules to current data frame
ess_ml_modules <- ess_ml %>%
  left_join(ise_gene_modules)

# remove module 0 as it is considered non-informative in the paper
ess_ml_modules <- ess_ml_modules %>% filter(cluster != 0)

# ============================================
# Fisher’s exact test and plot
# ============================================

# summarise the number of genes and number of 'peak' genes per module
modules_count <- ess_ml_modules %>%
  group_by(cluster) %>%
  summarise(
    total_genes = n(),
    peak_true   = sum(peak),
    .groups = "drop"
  )

# define background counts across all genes
background_total <- nrow(ess_ml)
background_peak  <- sum(ess_ml$peak)

# perform Fisher’s exact test for each module to test enrichment of 'peak' genes
fisher_results <- modules_count %>%
  rowwise() %>%
  mutate(
    a = peak_true,
    b = total_genes - a,
    c = background_peak - a,
    d = (background_total - total_genes) - c,
    test        = list(fisher.test(matrix(c(a,b,c,d), nrow = 2), alternative = "greater")),
    p_value     = test$p.value,
    odds_ratio  = as.numeric(test$estimate)
  ) %>%
  ungroup() %>%
  select(cluster, p_value, odds_ratio) %>%
  mutate(
    log2_OR     = log2(odds_ratio),
    neg_log10_p = -log10(p_value),
    is_significant = p_value < 0.05
  )

# Create the Volcano Plot showing module enrichment
ggplot(results, aes(x = log2_OR, y = neg_log10_p)) +
  geom_point(aes(color = p_value < 0.05), size = 6, shape = 19,  alpha = 0.8) +
  scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
                     labels = c("TRUE" = "P < 0.05", "FALSE" = "P >= 0.05"),
                     name = "Significance") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    x = "log2(Odds Ratio)",
    y = "-log10(p-value)",
    title = "Volcano Plot of modules Enrichment"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

