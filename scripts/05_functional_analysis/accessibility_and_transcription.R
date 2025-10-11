# ============================================
# Load R packages 
# ============================================
library(dplyr)
library(readr)

# ============================================
# Read inputs
# ============================================

# geneClusterCnts.tsv: 3 columns = ise, cluster, count (whitespace-separated)
gene_cluster <- read.table(
  "/scripts/05_functional_analysis/ise_geneClusterCnts.tsv",
  sep = "", header = FALSE,
  col.names = c("ise", "cluster", "count"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    # Normalise ISE ID: replace '-' with '_'
    ise = gsub("-", "_", ise),
    # Ensure integer count
    count = as.integer(count)
  )

# ISE gene names whitelist
ise_whitelist <- read.table(
  "/scripts/05_functional_analysis/ise_gene.txt",
  header = FALSE, col.names = c("ise")
)

# heacon-5 genes and transcripts names mapping
hc_ise <- read.csv(
  "data/heacon5_and_ise_mapping.csv"
)
hc_ise_pairs <- hc_ise %>%
  filter(Status != "Miss") %>%
  transmute(
    transcript_id = Query_ID,
    ise = gsub("^transcript:|(-\\d+$)", "", Subject_ID)
  )

# heacon-5 gene ↔ transcript mapping
hc_gene_transcript_mapping <- read.table(
  "/results/04_promoter_accessible_genes/heacon5_gene_and_transcript_mapping.txt",
  sep = "", header = FALSE,
  col.names = c("gene_id", "transcript_id")
)

# Promoter-accessible genes
pag <- read.table(
  "/results/04_promoter_accessible_genes/heacon5_promoter_accessible_genes.txt",
  header = FALSE, sep = "\t",
  col.names = c("gene_id", "transcript_id")
)

# ============================================
# Data frame praparation
# ============================================

# Complete cluster × ISE grid — fill zero for all missing combinations
all_clusters <- sort(unique(gene_cluster$cluster))
all_ise <- sort(unique(ise_whitelist$ise))

complete_combinations <- expand.grid(
  cluster = all_clusters,
  ise = all_ise,
  stringsAsFactors = FALSE
)

gene_cluster_completed <- complete_combinations %>%
  left_join(gene_cluster, by = c("cluster", "ise")) %>%
  mutate(count = ifelse(is.na(count), 0L, as.integer(count)))

# Attach mappings — ISE → heacon5 transcript_id, then heacon5 transcript_id → heacon5 gene
gene_cluster_completed <- gene_cluster_completed %>%
  left_join(hc_ise_pairs, by = "ise") %>%
  left_join(hc_gene_transcript_mapping, by = "transcript_id")

# ============================================
# Mark peak status — label using PAG gene_id set
# ============================================

# Gene IDs that are promoter-accessible (i.e., have peak)
pag_gene_set <- unique(pag$gene_id)

gene_cluster_completed <- gene_cluster_completed %>%
  mutate(
    peak = dplyr::case_when(
      is.na(gene_id) ~ NA,                  # no gene info → NA
      gene_id %in% pag_gene_set ~ TRUE,     # in PAG list → TRUE
      TRUE ~ FALSE                          # has gene but not in PAG → FALSE
    )
  )

# Keep rows that have a mapped gene (use this for per-gene stats)
gene_cluster_completed_keep <- gene_cluster_completed %>% filter(!is.na(gene_id))

# ============================================
# Summarise per gene — total count per gene; whether any peak exists
# ============================================

gene_summary <- gene_cluster_completed_keep %>%
  group_by(gene_id) %>%
  summarise(
    count = sum(count, na.rm = TRUE),
    peak  = any(peak == TRUE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(log_count = log10(count + 1))

# ============================================
# Overall stats & Wilcoxon test — TRUE vs FALSE (peak) global difference
# ============================================

overall_stats <- gene_summary %>%
  group_by(peak) %>%
  summarise(
    n            = n(),
    mean_count   = mean(count),
    median_count = median(count),
    sd_count     = sd(count),
    .groups = "drop"
  )

overall_wilcox <- wilcox.test(count ~ peak, data = gene_summary)
overall_stats$p_value   <- overall_wilcox$p.value
overall_stats$statistic <- unname(overall_wilcox$statistic)

# ============================================
# Per-cluster Wilcoxon — within each cluster: TRUE vs FALSE
# ============================================

per_cluster_wilcox <- gene_cluster_completed_keep %>%
  filter(!is.na(peak)) %>%
  group_by(cluster) %>%
  summarise(
    p_value      = tryCatch(wilcox.test(count ~ peak)$p.value, error = function(e) NA_real_),
    mean_TRUE    = mean(count[peak == TRUE],  na.rm = TRUE),
    mean_FALSE   = mean(count[peak == FALSE], na.rm = TRUE),
    median_TRUE  = median(count[peak == TRUE],  na.rm = TRUE),
    median_FALSE = median(count[peak == FALSE], na.rm = TRUE),
    n_TRUE       = sum(peak == TRUE,  na.rm = TRUE),
    n_FALSE      = sum(peak == FALSE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(cluster) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))  # optional: FDR correction