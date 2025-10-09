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
pag <- read.table("/results/04_promoter_accessible_genes/promoter_accessible_genes.txt",
                       header = FALSE,
                       sep = "\t",
                       stringsAsFactors = FALSE,
                       col.names = c("gene_id", "transcript_id"))

# ise and heacon-5 gene names mapping
hc_ise <- read.delim("/data/hc_ise_mapping.csv",
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
gene_map <- read.table("/data/hc_gene_transcript_mapping.txt",
                       sep = "",     
                       header = FALSE,
                       stringsAsFactors = FALSE,
                       col.names = c("gene_id", "transcript_id"))


top.E.gene <- read_excel(
  "/scripts/06_gene_essentiality_analysis/gene_essentiality_prediction_using_machine_learning.xlsx", col_names = c("ise", "rank")
)

tulio <- tulio %>%
  left_join(top.E.gene %>% dplyr::select(ise, rank),
            by = c("ise" = "ise"))

tulio <- tulio[!is.na(tulio$rank), ]
tulio <- tulio %>%
  mutate(peak = gene_id %in% gene_map$gene_id)

table(tulio$peak)

# Gene Essentiality prediction using Bingo
ess_bingo <- read.delim("/scripts/06_gene_essentiality_analysis/gene_essentiality_prediction_using_bingo.txt", header = TRUE, quote = "", comment.char = "")
ess_bingo <- ess_bingo %>%
  rename(ise = Transcript.ID) %>%
  mutate(ise = sub("-\\d+$", "", ise)) %>%
  mutate(rank = dense_rank(desc(Prediction.value)))

# add heacon-5 gene orthologs to bingo essentiality prediction
ess_bingo <- ess_bingo %>%
  left_join(hc_ise_pairs)

######## 导入 hc gene


gene_map.1 <- read.table("~/Desktop/research project/data/ref/hc_loner_protein_gene_ls.txt",
                         sep = "",           # 允许空格分隔
                         header = FALSE,     # 没有列名
                         stringsAsFactors = FALSE,

                                                  col.names = c("protein", "gene_id"))
jiani <- jiani %>%
  left_join(gene_map)

jiani <- jiani %>%
  left_join(gene_map.1)


#### peak

gene.map <- read.table("~/Desktop/research project/data/tss_summit_gene_map.2.txt",
                       header = FALSE,
                       sep = "\t",
                       stringsAsFactors = FALSE,
                       col.names = c("gene_id", "transcript_id"))

jiani <- jiani %>%
  mutate(peak = ifelse(is.na(gene_id), NA, gene_id %in% gene.map$gene_id))


write.table(jiani, file = "Desktop/research project/jiani/jiani.tsv", sep = "\t", quote = FALSE, row.names = FALSE)





# # filter group 0 out
# tulio_filtered <- tulio %>% filter(cluster != 0)
# 
# # find enriched gene from scrna.R
# tulio_filtered <- tulio_filtered %>%
#   dplyr::mutate(`isfunction` = protein %in% final_df$gene)
# 
# # find is top1000 or not
# tulio_filtered %>%
#   group_by(cluster, isTop) %>%
#   summarise(count = n(), .groups = "drop")
# 
# 
# 
# 
# 
# 
# # 创建列联表
# table_cluster_isTop <- table(tulio_filtered$cluster, tulio_filtered$isTop)
# 
# # 执行卡方检验
# chisq.test(table_cluster_isTop)



# # 先确认 rank 是数值型
# str(tulio_filtered$rank)
# 
# # 进行单因素方差分析（ANOVA）
# anova_result <- aov(rank ~ as.factor(cluster), data = tulio_filtered)
# summary(anova_result)


# ggplot(tulio_filtered, aes(x = factor(cluster), y = rank)) +  # Use 'rank' for the y-value
#   geom_boxplot(fill = "lightblue") +
#   scale_y_reverse() +  # Reverses the y-axis: rank 1 will be higher up
#   xlab("Cluster") +    # X-axis label
#   ylab("essentiality") + # Y-axis title displayed as "essentiality"
#   ggtitle("Distribution of essentiality (Rank) by Cluster") + # Updated plot title
#   theme(
#     axis.text.y = element_blank(),  # Remove y-axis tick labels (numerical values)
#     axis.ticks.y = element_blank(), # Remove y-axis tick marks
#     axis.line.y = element_blank(),  # Remove the y-axis line
#     # If you also want to remove the y-axis title "essentiality", add:
#     # axis.title.y = element_blank(),
#     panel.background = element_rect(fill = "white"), # Optional: set background to white
#     panel.grid.major = element_line(colour = "grey97"), # Optional: light grid lines
#     panel.grid.minor = element_blank() # Optional: remove minor grid lines
#   )
# 
# 



library(dplyr)
library(ggplot2)
library(tidyr)

summary_df <- tulio_filtered %>%
  group_by(cluster) %>%
  summarise(
    total = n(),
    peak_ratio = mean(peak, na.rm = TRUE),
    function_ratio = mean(isfunction, na.rm = TRUE),
    top_ratio = mean(isTop, na.rm = TRUE),
    .groups = "drop"
  )

# 转换成长格式以便 ggplot 绘图
long_df <- summary_df %>%
  pivot_longer(cols = c(peak_ratio, function_ratio, top_ratio),
               names_to = "type", values_to = "ratio")


# peak vs isTop
fisher.test(table(tulio_filtered$peak, tulio_filtered$isTop))
# isfunction vs isTop
fisher.test(table(tulio_filtered$isfunction, tulio_filtered$isTop))






heat <- tulio_filtered[, c("cluster", "peak", "isfunction")]

head(heat)


# 汇总统计
cluster_df <- tulio_filtered %>%
  group_by(cluster) %>%
  summarise(
    total_genes = n(),
    peak_true = sum(peak == TRUE)
  )

sum(cluster_df$total_genes)
sum(cluster_df$peak_true)

# 
# cluster_tulio_summary <- tulio %>%
#   group_by(cluster) %>%
#   summarise(
#     total_genes = n(),
#     peak_true = sum(peak == TRUE)
#   )
# 
# sum(cluster_tulio_summary$total_genes)
# sum(cluster_tulio_summary$peak_true)

# 假设你的数据叫 cluster_df，包含 cluster, total_genes, peak_true
background_total <- 14825
background_peak <- 3253

results <- lapply(1:nrow(cluster_df), function(i) {
  a <- cluster_df$peak_true[i]
  b <- cluster_df$total_genes[i] - a
  c <- background_peak - a
  d <- (background_total - cluster_df$total_genes[i]) - c
  
  mat <- matrix(c(a, b, c, d), nrow = 2)
  fisher <- fisher.test(mat, alternative = "greater")  # 检测是否 enrichment
  
  data.frame(
    cluster = cluster_df$cluster[i],
    p_value = fisher$p.value,
    odds_ratio = fisher$estimate
  )
}) %>% dplyr::bind_rows()


results <- results %>%
  mutate(
    `-log10(p_value)` = -log10(p_value),
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Enrichment = ifelse(odds_ratio > 1, "Enriched", "Not enriched")
  )



results <- results %>%
  mutate(
    log2_OR = log2(odds_ratio),
    neg_log10_p = -log10(p_value),
    is_significant = p_value < 0.05 # Create a boolean column for significance
  )

# Define the significance threshold for the p-value
p_threshold <- 0.05


# fisher-test for genes' essentiality and accessibility
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
    title = "Volcano Plot of Cluster Enrichment"
  ) +
  theme_minimal() + # Still uses a clean base theme
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    # --- Add this line to draw the axis lines ---
    axis.line = element_line(color = "black", linewidth = 0.5)
    # ---------------------------------------------
  )








str(tulio)

df_ranked <- tulio_filtered %>% filter(!is.na(rank))









ks.test(rank ~ peak, data = df_ranked)


ggplot(df_ranked, aes(x = rank, color = peak)) +
  # 使用 stat_ecdf() 绘制经验累积分布曲线
  stat_ecdf(geom = "step", size = 1) +
  
  # 手动设置颜色, 与您之前的例子保持一致
  scale_color_manual(
    name = "Accessibility", # 设置图例标题
    values = c("FALSE" = "gray50", "TRUE" = "red3"),
    labels = c("FALSE" = "Inaccessible", "TRUE" = "Accessible")
  ) +
  
  # 设置坐标轴和标题
  labs(
    x = "Essentiality Rank",
    y = "Cumulative Probability",
    title = "ECDF of Essentiality Rank by Accessibility",
    subtitle = "Visualizing distributions for the Kolmogorov-Smirnov test"
  ) +
  
  # 使用简洁主题
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom", # 将图例放在底部
    legend.title = element_text(face = "bold"),
    panel.grid.minor.x = element_blank(),  # 移除 y 轴次线
  )









ks_res <- ks.test(df_ranked$rank[df_ranked$peak == TRUE],
                  df_ranked$rank[df_ranked$peak == FALSE])










ks_result <- ks.test(rank ~ peak, data = df_ranked)




# 假设之前的代码（数据创建、KS检验）已经运行

ggplot(df_ranked, aes(x = peak, y = rank, fill = peak)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, color = "black", outlier.size = 0.5, alpha = 0.7) +
  scale_y_reverse() +
  scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "red3")) +
  scale_x_discrete(labels = c("FALSE" = "Inaccessible", "TRUE" = "Accessible")) +
  # 移除坐标轴留白
  scale_y_continuous(expand = expansion(mult = c(0.05, 0))) + 
  
  labs(
    x = NULL,
    y = "Essentiality Rank",
    title = "Essentiality Rank vs Accessibility",
    subtitle = paste0("KS test: D = ", ks_d, ", ", ks_p)
  ) +
  theme_minimal(base_size = 14) +
  
  # --- 这是修改的部分 ---
  theme(
    panel.grid.major.y = element_blank(),  # 移除 y 轴主线
    panel.grid.minor.y = element_blank(),  # 移除 y 轴次线
    panel.grid.minor.x = element_blank(),  # 移除 y 轴次线
    legend.position = "none",
    axis.text.x = element_text(size = 12), 
  ) +
  coord_flip()













# Bar plot of -log10(p_value)
ggplot(results, aes(x = reorder(cluster, `-log10(p_value)`), y = `-log10(p_value)`, fill = Enrichment)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(label = Significance), vjust = -0.5, size = 5) + # Add significance stars
  scale_fill_manual(values = c("Enriched" = "skyblue", "Not enriched" = "lightcoral")) +
  labs(
    title = "Enrichment Analysis of Clusters",
    x = "Cluster",
    y = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Bar plot of Odds Ratio
ggplot(results, aes(x = reorder(cluster, odds_ratio), y = odds_ratio, fill = Enrichment)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") + # Reference line for odds ratio = 1
  scale_fill_manual(values = c("Enriched" = "skyblue", "Not enriched" = "lightcoral")) +
  labs(
    title = "Odds Ratio for Cluster Enrichment",
    x = "Cluster",
    y = "Odds Ratio"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



















































