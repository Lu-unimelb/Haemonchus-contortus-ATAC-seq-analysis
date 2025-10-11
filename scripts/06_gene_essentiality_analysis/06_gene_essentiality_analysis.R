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
# promoter accessible gene
pag <- read.table("/results/04_promoter_accessible_genes/heacon5_promoter_accessible_genes.txt",
                       header = FALSE,
                       sep = "\t",
                       stringsAsFactors = FALSE,
                       col.names = c("gene_id", "transcript_id"))

# ise and heacon-5 gene names mapping
hc_ise <- read.delim("/data/heacon5_and_ise_mapping.csv",
                     header = TRUE, stringsAsFactors = FALSE)
hc_ise_pairs <- hc_ise %>%
  filter(Status != "Miss") %>%
  transmute(
    hc  = Query_ID,
    ise = gsub("^transcript:|(-\\d+$)", "", Subject_ID)
  )
colnames(hc_ise_pairs) <- c("transcript_id", "ise")
hc_ise_pairs <- hc_ise_pairs %>%
  filter(!is.na(transcript_id) & transcript_id != "null")

# heacon-5 genes and transcripts names mapping
hc_gene_transcript_mapping <- read.table("/data/heacon5_gene_and_transcript_mapping.txt",
                       sep = "",     
                       header = FALSE,
                       stringsAsFactors = FALSE,
                       col.names = c("gene_id", "transcript_id"))

# heacon-5 longest protein isoform and gene mapping names mapping
hc_protein_gene_mapping <- read.table("/data/heacon5_longest_protein_isoform_and_gene_mapping.txt",
                         sep = "",
                         header = FALSE,
                         stringsAsFactors = FALSE,
                         col.names = c("protein", "gene_id"))

# ============================================
# Accessibility and Essentiality (BINGO) data frame construction
# ============================================

# Gene Essentiality prediction using Bingo
ess_bingo <- read.delim("/scripts/06_gene_essentiality_analysis/gene_essentiality_prediction_using_bingo.txt", header = TRUE, quote = "", comment.char = "")
ess_bingo <- ess_bingo %>%
  rename(ise = Transcript.ID) %>%
  mutate(ise = sub("-\\d+$", "", ise)) %>%
  mutate(rank = dense_rank(desc(Prediction.value)))

# add heacon-5 gene orthologs to bingo essentiality prediction data frame
ess_bingo <- ess_bingo %>%
  left_join(hc_ise_pairs)

# add heacon-5 transcripts to data frame
ess_bingo <- ess_bingo %>%
  left_join(hc_gene_transcript_mapping)

# add heacon-5 protein to data frame
ess_bingo <- ess_bingo %>%
  left_join(hc_protein_gene_mapping)

# mark promoter accessible genes
ess_bingo <- ess_bingo %>%
  mutate(peak = ifelse(is.na(gene_id), NA, gene_id %in% pag$gene_id))

# label essentiality rank to each gene
ess_bingo$rank <- seq_len(nrow(ess_bingo))

# remove genes without heacon5 orthologs
ess_bingo <- ess_bingo %>% filter(!is.na(peak))

# ============================================
# KS Test of Accessibility and Essentiality (BINGO) and plot
# ============================================

# compare accessible vs. inaccessible genes' essentiality ranks using the Kolmogorov–Smirnov test
ks.test(rank ~ peak, data = ess_bingo)

# plot the empirical cumulative distribution (ECDF) of essentiality ranks
ggplot(ess_bingo, aes(x = rank, color = peak)) +
  stat_ecdf(geom = "step", size = 1) +

  scale_color_manual(
    name = "Accessibility",
    values = c("FALSE" = "gray50", "TRUE" = "red3"),
    labels = c("FALSE" = "Inaccessible", "TRUE" = "Accessible")
  ) +

  labs(
    x = "Essentiality Rank",
    y = "Cumulative Probability",
    title = "ECDF of Essentiality Rank by Accessibility",
    subtitle = "Visualizing distributions for the Kolmogorov-Smirnov test"
  ) +

  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor.x = element_blank(),
  )

# Visualise the distribution of essentiality ranks between accessible and inaccessible genes using a violin plot
ggplot(ess_bingo, aes(x = peak, y = rank, fill = peak)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, color = "black", outlier.size = 0.5, alpha = 0.7) +
  scale_y_reverse() +
  scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "red3")) +
  scale_x_discrete(labels = c("FALSE" = "Inaccessible", "TRUE" = "Accessible")) +

  scale_y_continuous(expand = expansion(mult = c(0.05, 0))) + 
  
  labs(
    x = NULL,
    y = "Essentiality Rank",
    title = "Essentiality Rank vs Accessibility"
  ) +
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 12), 
  ) +
  coord_flip()

# ============================================
# Accessibility and Essentiality (ML) data frame construction
# ============================================

ess_ml <- read_excel(
  "gene_essentiality_prediction_using_machine_learning.xlsx",
  sheet = 1, col_names = c("ise", "rank")
)

# add heacon-5 gene orthologs to ml essentiality prediction data frame
ess_ml <- ess_ml %>%
  left_join(hc_ise_pairs)

# add heacon-5 transcripts to data frame
ess_ml <- ess_ml %>%
  left_join(hc_gene_transcript_mapping)

# add heacon-5 protein to data frame
ess_ml <- ess_ml %>%
  left_join(hc_protein_gene_mapping)

# mark promoter accessible genes
ess_ml <- ess_ml %>%
  mutate(peak = ifelse(is.na(gene_id), NA, gene_id %in% pag$gene_id))

# remove genes without heacon5 orthologs
ess_ml <- ess_ml %>% filter(!is.na(peak))

# save data frame for next step
write.table(ess_ml, file = "/results/06_gene_essentiality_analysis/ess_ml.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================
# KS Test of Accessibility and Essentiality (ML) and plot
# ============================================

# compare accessible vs. inaccessible genes' essentiality ranks using the Kolmogorov–Smirnov test
ks.test(rank ~ peak, data = ess_ml)

# plot the empirical cumulative distribution (ECDF) of essentiality ranks
ggplot(ess_ml, aes(x = rank, color = peak)) +
  stat_ecdf(geom = "step", size = 1) +
  
  scale_color_manual(
    name = "Accessibility",
    values = c("FALSE" = "gray50", "TRUE" = "red3"),
    labels = c("FALSE" = "Inaccessible", "TRUE" = "Accessible")
  ) +
  
  labs(
    x = "Essentiality Rank",
    y = "Cumulative Probability",
    title = "ECDF of Essentiality Rank by Accessibility",
    subtitle = "Visualizing distributions for the Kolmogorov-Smirnov test"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor.x = element_blank(),
  )

# Visualise the distribution of essentiality ranks between accessible and inaccessible genes using a violin plot
ggplot(ess_ml, aes(x = peak, y = rank, fill = peak)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, color = "black", outlier.size = 0.5, alpha = 0.7) +
  scale_y_reverse() +
  scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "red3")) +
  scale_x_discrete(labels = c("FALSE" = "Inaccessible", "TRUE" = "Accessible")) +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0))) + 
  
  labs(
    x = NULL,
    y = "Essentiality Rank",
    title = "Essentiality Rank vs Accessibility"
  ) +
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 12), 
  ) +
  coord_flip()



















