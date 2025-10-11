# ============================================
# Load R packages 
# ============================================
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(clusterProfiler)
library(GO.db)
library(enrichplot)
library(KEGGREST)
library(forcats)

# ============================================
# File preparation
# ============================================
# background genes
universe <- read.delim("/data/heacon5_longest_protein_isoform.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1]
# Go term to gene 
go_annotations <- read.table("/scripts/05_functional_analysis/enrichment_inputs/GO_annotation.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(go_annotations) <- c("gene", "GOs")
# Kegg orthologs to gene
ko2gene <- read.table("/scripts/05_functional_analysis/enrichment_inputs/ko2gene.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(ko2gene) <- c("gene", "ko")
# Kegg orthologs description
ko <- read.delim("/scripts/05_functional_analysis/enrichment_inputs/ko2name.tsv",
                 header     = FALSE,
                 sep        = "\t",
                 col.names  = c("ko","Description"),
                 colClasses = c("character","character"),
                 stringsAsFactors = FALSE)
# Kegg pathway to gene
path2gene <- read.table("/scripts/05_functional_analysis/enrichment_inputs/path2gene.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(path2gene) <- c("gene", "path")
# Kegg pathway description 
kpath <- read.delim("/scripts/05_functional_analysis/enrichment_inputs/path2name.tsv", 
                    header     = FALSE,
                    sep        = "\t",
                    col.names  = c("path","Description"),
                    colClasses = c("character","character"),
                    stringsAsFactors = FALSE)
# Go Terms
go2name <- Term(GOTERM)
go2name=as.data.frame(go2name)
go2name <- rownames_to_column(go2name,var = 'GO_id') %>% 
  mutate(GO_id = str_sub(GO_id,4,10) )

# Data.frame cleaning
go2gene = go_annotations %>% 
  dplyr::select(1,2) %>% 
  filter(GOs != "-") %>% 
  separate_rows(GOs,sep = ",",convert = F) %>% 
  mutate(GOs = substr(GOs,4,10)) %>% 
  unique() %>% 
  dplyr::select(2,1)
ko2gene <- ko2gene %>%
  filter(ko != "-") %>%
  separate_rows(ko, sep = ",",convert = F) %>%
  mutate(ko = str_remove(ko, "ko:")) %>%
  unique() %>% dplyr::select(2,1)
# Pathways data.frame cleaning
path2gene <- path2gene %>%
  filter(path != "-") %>%
  mutate(path = str_sub(path, -5)) %>%
  unique() %>% dplyr::select(2,1)
# Pathway description
path_id = c("04010", "04510", "04015", "04151", "04360", "04024", "04014", "04020", "04022", "04072", "04012", "04211", "04810", "04390", "04724", "04725")
path_d = c("MAPK signaling pathway", "Focal adhesion", "Rap1 signaling pathway", "PI3K-Akt signaling pathway", "Axon guidance", "cAMP signaling pathway", "Ras signaling pathway", "Calcium signaling pathway", "cGMP-PKG signaling pathway", "Phospholipase D signaling pathway", "ErbB signaling pathway", "	Longevity regulating pathway", "Regulation of actin cytoskeleton", "Hippo signaling pathway", "Glutamatergic synapse", "Cholinergic synapse")
path_new <- data.frame(
  path = path_id,
  Description = path_d,
  stringsAsFactors = FALSE
)
kpath <- rbind(kpath, path_new) %>% arrange(path)

kegg_term <- read_tsv(
  "/scripts/05_functional_analysis/enrichment_inputs/kegg_term.txt",
  col_types = cols(
    .default = col_character()
  )
)

# promoter accessible gene
pag <- read.table("/results/04_promoter_accessible_genes/heacon5_promoter_accessible_genes.txt", 
                       header = FALSE, 
                       sep = "\t", 
                       stringsAsFactors = FALSE,
                       col.names = c("gene_id", "transcript_id"))

pag_df <- as.data.frame(table(pag$gene_id))
colnames(pag_df) <- c("gene_id", "transcript_count")

protein_gene_mapping <- read.table(
  "/data/heacon5_longest_protein_isoform_and_gene_mapping.txt",
  header = FALSE, stringsAsFactors = FALSE
)
colnames(protein_gene_mapping) <- c("protein", "gene")

pag_df <- pag_df %>%
  left_join(protein_gene_mapping %>% dplyr::select(gene, protein),
            by = c("gene_id" = "gene"))

key_gene <- pag_df$protein %>% as.character()

# ============================================
# Enrichment and plot results function
# ============================================

plot_enrichment <- function(key_gene, universe, 
                                go2gene, go2name,
                                path2gene, kpath,
                                kegg_term,
                                title = "Enriched Terms") {
  
  # --- GO enrichment ---
  go_top10_x <- NULL
  go_result <- tryCatch({
    enricher(key_gene, universe = universe, TERM2GENE = go2gene, TERM2NAME = go2name,
             pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  }, error = function(e) NULL)
  
  if (!is.null(go_result) && !is.null(go_result@result) && nrow(go_result@result) > 0) {
    go_df <- as.data.frame(go_result@result) %>% 
      filter(Description != "NA", p.adjust < 0.05)
    
    go_ont <- go2ont(str_c("GO:", go_df$ID)) %>% 
      mutate(go_id = str_sub(go_id, 4, 10))
    
    go_df <- left_join(go_df, go_ont, by = c("ID" = "go_id")) %>% 
      filter(Ontology == "BP") %>% 
      arrange(p.adjust)
    
    go_top10 <- go_df %>%
      group_by(Ontology) %>%
      slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>%
      ungroup()
    
    go_top10$Short_Description <- sapply(go_top10$Description, function(x) {
      words <- strsplit(as.character(x), split = " ")[[1]]
      paste(words[1:min(10, length(words))], collapse = " ")
    })
    
    go_top10 <- go_top10 %>%
      mutate(log10p = -log10(p.adjust)) %>%
      group_by(Ontology) %>%
      arrange(desc(log10p), .by_group = TRUE) %>%
      ungroup() %>%
      mutate(type_order = factor(Short_Description, levels = rev(unique(Short_Description))))
    
    go_top10_x <- go_top10 %>% dplyr::select(type_order, log10p, Ontology)
  } else {
    message("No GO enrichment results.")
  }
  
  
  # --- Pathway enrichment ---
  path_top10_x <- NULL
  path_result <- tryCatch({
    enricher(key_gene, universe = universe, TERM2GENE = path2gene, TERM2NAME = kpath,
             pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  }, error = function(e) NULL)
  
  if (!is.null(path_result) && !is.null(path_result@result) && nrow(path_result@result) > 0) {
    path_df <- as.data.frame(path_result@result) %>% 
      filter(Description != "NA", p.adjust < 0.05)
    
    path_df <- path_df[!grepl("^\\d{5}$", path_df$Description), ]
    
    if (nrow(path_df) > 0) {
      path_top10 <- path_df %>% arrange(p.adjust) %>% slice_min(order_by = p.adjust, n = 30)

      path_top10$Short_Description <- path_top10$Description
      
      path_top10 <- path_top10 %>%
        mutate(log10p = -log10(p.adjust))

      A_level$row <- seq_len(nrow(A_level))
      C_rows <- which(A_level$Level == "C")
      
      C2A_map <- lapply(C_rows, function(i) {
        c_code <- A_level$Code[i]
        a_row <- tail(which(A_level$Level[1:i] == "A"), 1)
        list(Code = c_code, A_Description = A_level$Description[a_row])
      }) %>% bind_rows()

      path_top10 <- left_join(path_top10, C2A_map, by = c("ID" = "Code")) %>%
        mutate(Ontology = A_Description)

      path_top10_x <- path_top10 %>%
        dplyr::select(Short_Description, log10p, Ontology) %>%
        dplyr::rename(type_order = Short_Description)
    } else {
      message("No KEGG terms after filtering.")
    }
  } else {
    message("No KEGG Pathway enrichment results.")
  }
  
  # --- Combine and plot ---
  combined_data <- bind_rows(go_top10_x, path_top10_x)
  
  if (nrow(combined_data) == 0) {
    message("No enrichment results to plot.")
    return(NULL)
  }
  
  # --- Data Preparation ---
  # Ensure Ontology is a factor
  combined_data$Ontology <- factor(combined_data$Ontology)
  
  # Recode "BP" to "Biological Process" if "BP" is a level in Ontology
  if ("BP" %in% levels(combined_data$Ontology)) {
    combined_data$Ontology <- forcats::fct_recode(combined_data$Ontology, 
                                                  "Biological Process" = "BP")
  }
  
  # Define custom lighter colors - ensure all Ontology levels in data have a color
  # Get unique ontology levels from the data AFTER potential recoding
  present_ontologies <- levels(combined_data$Ontology)
  
  defined_colors <- c(
    "Biological Process" = "#fddede", # Light red/pink
    "Cellular Processes" = "#e3b7b7", # Desaturated red
    "Environmental Information Processing" = "#b3e5f7", # Light blue
    "Genetic Information Processing" = "#d8c3eb", # Light purple
    "Metabolism" = "#b3c8f4", # Light desaturated blue/purple
    "Organismal Systems" = "#f5d1f2", # Light pink/magenta
    "Pathway" = "#d3d3d3" # Default light grey for "Pathway" or other categories
  )
  
  # Filter custom_colors to only those present, and add defaults for any missing
  custom_colors_final <- defined_colors[names(defined_colors) %in% present_ontologies]
  missing_ontologies <- present_ontologies[!present_ontologies %in% names(custom_colors_final)]
  if(length(missing_ontologies) > 0) {
    # Create default colors for missing ontologies (e.g., shades of grey or other distinct colors)
    # For simplicity, using a generic palette here. You might want a more sophisticated approach.
    default_palette <- RColorBrewer::brewer.pal(n = max(3, length(missing_ontologies)), name = "Pastel2")
    if (length(missing_ontologies) > length(default_palette)) {
      default_palette <- rep(default_palette, length.out = length(missing_ontologies))
    }
    names(default_palette) <- missing_ontologies
    custom_colors_final <- c(custom_colors_final, default_palette[1:length(missing_ontologies)])
  }
  
  
  # Sort terms within each Ontology group by log10p
  # fct_reorder makes terms with higher log10p appear higher on the y-axis
  combined_data <- combined_data %>%
    dplyr::filter(!is.na(log10p) & !is.na(type_order) & !is.na(Ontology)) %>% # Ensure no NAs that break grouping/ordering
    dplyr::group_by(Ontology) %>%
    # Ensure type_order is character before fct_reorder if it became factor globally
    dplyr::mutate(type_order = as.character(type_order), 
                  type_order = forcats::fct_reorder(type_order, log10p)) %>%
    dplyr::ungroup()
  
  if (nrow(combined_data) == 0) {
    message("No valid data remaining after preparation for plotting.")
    return(NULL)
  }
  
  # --- Create the Plot ---
  p <- ggplot2::ggplot(combined_data, ggplot2::aes(x = log10p, y = type_order, fill = Ontology)) +
    ggplot2::geom_col(width = 0.8) + 
    ggplot2::geom_text(ggplot2::aes(x = 0.1, label = type_order),
                       hjust = 0, 
                       color = "black",
                       size = 3.2) +
    ggplot2::scale_fill_manual(values = custom_colors_final) +
    ggplot2::xlab(expression(-log[10]("adjusted p-value"))) +
    ggplot2::ylab(NULL) + 
    ggplot2::labs(title = title) +
    
    # MODIFIED PART: Using ggforce::facet_col for consistent bar height across facets
    # This requires the ggforce package: install.packages("ggforce")
    # It makes panel heights proportional to the number of items in their y-scale,
    # so geom_col(width=0.8) results in bars of the same absolute thickness.
    ggforce::facet_col(ggplot2::vars(Ontology), scales = "free_y", space = "free", strip.position = "top") +
    # ALTERNATIVE (if ggforce is not available, uncomment below and comment out facet_col):
    # This will make bar heights consistent but might add empty space in facets with fewer items.
    # ggplot2::facet_wrap(~Ontology, ncol = 1, strip.position = "top") + # Removed scales = "free_y"
    
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 11, margin = ggplot2::margin(t = 10)),
      axis.text.x = ggplot2::element_text(size = 10, color = "black"), 
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.ticks.x = ggplot2::element_line(color = "black", linewidth = 0.5),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14, margin = ggplot2::margin(b = 15)),
      legend.position = "none",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "grey95", linewidth = 0.5),
      panel.grid.minor.x = ggplot2::element_blank(), 
      strip.background = ggplot2::element_blank(), 
      strip.text = ggplot2::element_text( 
        hjust = 0, 
        face = "bold",
        size = 10,
        margin = ggplot2::margin(t = 1, b = 3) 
      ),
      panel.spacing.y = ggplot2::unit(0, "lines") 
    ) +
    ggplot2::coord_cartesian(clip = "off") 
  
  
  return(p)
}

# ============================================
# Running function
# ============================================

plot_enrichment_all(key_gene = key_gene,
                    universe = universe,
                    go2gene = go2gene,
                    go2name = go2name,
                    path2gene = path2gene,
                    kpath = kpath,
                    kegg_term = kegg_term,
                    title = "Enrichment")
