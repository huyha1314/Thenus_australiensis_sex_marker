#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(VennDiagram)
  library(patchwork)
  library(grid)
  library(stringr)
})

# ==============================================================================
# 0. COMMAND LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path to master_annotation_summary.tsv", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="Figure2_Combined_Final", 
              help="Prefix for output files", metavar="character")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ==============================================================================
# 1. GENERATE PANEL A: VENN DIAGRAM
# ==============================================================================
data <- read.delim(opt$input, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")

get_genes <- function(df, col) {
  df %>% filter(!!sym(col) != "-" & !!sym(col) != "") %>% pull(Gene_ID) %>% unique()
}

# The lists are mapped to ellipses 1, 2, 3, and 4 in this exact order
gene_lists <- list(
  KEGG = get_genes(data, "KEGG_ko"),                 # Ellipse 1: Bottom-Left
  NR = get_genes(col = "NR_Hit_ID", df = data),      # Ellipse 2: Bottom-Right
  SwissProt = get_genes(data, "SwissProt_Hit_ID"),   # Ellipse 3: Top-Left
  InterPro = get_genes(data, "InterPro_Acc")         # Ellipse 4: Top-Right
)

cb_palette <- c("#56B4E9", "#E69F00", "#009E73", "#CC79A7")

p_venn <- grid.grabExpr({
  grid.draw(venn.diagram(
    x = gene_lists, filename = NULL,
    fill = cb_palette, 
    alpha = 0.2,           
    col = "black", lwd = 2,
    cex = 2.6, fontface = "bold", fontfamily = "sans",
    cat.cex = 3.0, cat.fontface = "bold", cat.fontfamily = "sans",
    
    # THE FIX: Properly mapping the angles to the corresponding ellipses
    # Order: Bottom-Left (-135), Bottom-Right (135), Top-Left (-45), Top-Right (45)
    cat.pos = c(-135, 135, -45, 45), 
    
    # Brought the distances in slightly so they stay neatly inside the bounding box
    cat.dist = c(0.22, 0.22, 0.22, 0.22),
    margin = 0.25
  ))
})

# ==============================================================================
# 2. GENERATE PANEL B: COG BARPLOT
# ==============================================================================
cog_raw <- data %>% filter(COG_category != "-" & COG_category != "") %>% pull(COG_category)
cog_counts <- data.frame(Category = unlist(strsplit(cog_raw, ""))) %>%
  group_by(Category) %>% summarise(Count = n(), .groups = 'drop')

cog_key <- read.delim("https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.fun.tab", 
                      header = FALSE, sep = "\t", comment.char = "#")
colnames(cog_key) <- c("Category", "FunctionalGroup", "ColorHex", "Description")

plot_data <- merge(cog_counts, cog_key, by ='Category') %>%
  mutate(Group = case_when(
    Category %in% c("J", "A", "K", "L", "B") ~ "Information storage and processing",
    Category %in% c("D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O") ~ "Cellular processes and signaling",
    Category %in% c("C", "G", "E", "F", "H", "I", "P", "Q") ~ "Metabolism",
    Category %in% c("R", "S") ~ "Poorly characterized",
    TRUE ~ "Other"
  )) %>%
  arrange(Group, Count) %>%
  mutate(Category = factor(Category, levels = Category),
         key_df = str_wrap(paste0("[", Category, "] ", Description), width = 28))

p_cog_bar <- ggplot(plot_data, aes(x = Category, y = Count, fill = Group)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.8) +
  geom_text(aes(label = Count), vjust = -0.5, size = 9, fontface="bold") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 42),
        axis.title = element_text(size = 34, face = "bold"),
        axis.text = element_text(size = 29, color = "black"),
        legend.position = "top", 
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 29, face = "bold"))

# ==============================================================================
# 3. GENERATE PANEL C: 4-COLUMN LEGEND
# ==============================================================================
key_data <- plot_data %>% arrange(Group, Category) %>% pull(key_df)

num_cols <- 4
num_rows <- ceiling(length(key_data) / num_cols)
key_data_padded <- c(key_data, rep("", (num_rows * num_cols) - length(key_data)))

col_matrix <- matrix(key_data_padded, nrow = num_rows, ncol = num_cols, byrow = FALSE)
key_df_fmt <- as.data.frame(col_matrix)

y_coords <- seq(1, by = 4.5, length.out = num_rows)

p_key <- ggplot(key_df_fmt) + theme_void() +
  geom_text(aes(x = 0.00, y = y_coords, label = V1), hjust = 0, vjust = 1, size = 8.5, lineheight = 0.85) +
  geom_text(aes(x = 0.25, y = y_coords, label = V2), hjust = 0, vjust = 1, size = 8.5, lineheight = 0.85) +
  geom_text(aes(x = 0.50, y = y_coords, label = V3), hjust = 0, vjust = 1, size = 8.5, lineheight = 0.85) +
  geom_text(aes(x = 0.75, y = y_coords, label = V4), hjust = 0, vjust = 1, size = 8.5, lineheight = 0.85) +
  scale_y_reverse(limits = c(max(y_coords) + 5, 0)) + xlim(0, 1)

# ==============================================================================
# 4. COMBINE AND SAVE
# ==============================================================================
final_figure <- (wrap_elements(p_venn) / p_cog_bar / p_key) + 
  plot_layout(heights = c(1.5, 1.2, 1.2)) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') & 
  theme(plot.tag = element_text(size = 45, face = "bold"))

ggsave(paste0(opt$out, ".png"), plot = final_figure, width = 26, height = 32, dpi = 300, bg = "white")
ggsave(paste0(opt$out, ".pdf"), plot = final_figure, width = 26, height = 32, device = cairo_pdf)