#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# 0. COMMAND LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  make_option(c("-m", "--matrix"), type="character", default="cenpe-expr.matrix", 
              help="Path to expression matrix file [default= %default]", metavar="character"),
  make_option(c("-s", "--meta"), type="character", default="sample.tsv", 
              help="Path to metadata file [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="Figure_Combined_CENPE_Profile", 
              help="Prefix for output files [default= %default]", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Verify inputs exist before proceeding
if (!file.exists(opt$matrix)) stop(paste("Error: Matrix file not found at", opt$matrix))
if (!file.exists(opt$meta)) stop(paste("Error: Metadata file not found at", opt$meta))

# ==============================================================================
# 1. SETUP: COLOR PALETTE & THEME
# ==============================================================================
# Okabe-Ito Colorblind-Safe Palette
sex_colors <- c("Male" = "#56B4E9", "Female" = "#D55E00")

# Colorblind-safe replacements for the Red/Green protein length bars
len_colors <- c("Ancestral" = "#2A4B7C", "NeoY" = "#009E73")

# Define a consistent publication-ready Theme
publication_theme <- theme_classic() +
  theme(
    plot.title = element_text(size = 26, face = "bold"),
    plot.subtitle = element_text(size = 20, color = "grey30", margin = margin(b=15)),
    axis.title.x = element_text(size = 22, face = "bold", margin = margin(t=15)),
    axis.title.y = element_text(size = 22, face = "bold", margin = margin(r=15)),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1, face = "bold", color="black"),
    axis.text.y = element_text(size = 18, color="black"),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    legend.position = "top",
    panel.grid.major.y = element_line(color="grey90", linetype="dashed")
  )

# ==============================================================================
# 2. DEFINE PLOTTING FUNCTIONS
# ==============================================================================

# FUNCTION A: Expression Bar Plot
generate_tmm_plot <- function(gene_name, summary_df) {
  plot_data <- summary_df %>% filter(gene_id == gene_name)
  
  ggplot(plot_data, aes(x = tissue, y = mean_tmm, fill = sex)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth=0.8) +
    geom_errorbar(aes(ymin = mean_tmm - se_tmm, ymax = mean_tmm + se_tmm),
                  position = position_dodge(width = 0.8), width = 0.25, linewidth=0.8) +
    scale_fill_manual(values = sex_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
    labs(title = paste("Expression Profile:", gene_name),
         subtitle = "Mean TMM ± Standard Error",
         y = "Mean TMM Expression",
         x = "Tissue",
         fill = "Sex") +
    publication_theme
}

# FUNCTION B: Protein Length Horizontal Bar
generate_length_plot <- function(transcript_id, prot_length, bar_color) {
  df <- data.frame(id = transcript_id, length = prot_length)
  
  ggplot(df, aes(x = length, y = id)) +
    geom_col(fill = bar_color, color = "black", width = 0.3, linewidth=1) +
    # Lock the X-axis from 0 to 3200 so both plots align perfectly visually
    scale_x_continuous(limits = c(0, 3200), breaks = seq(0, 3000, by=500), expand = c(0,0)) +
    labs(x = "Protein Length (aa)", y = NULL) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 20, face = "bold", margin = margin(t=10)),
      axis.text.x = element_text(size = 18, color = "black"),
      # Rotate the Y-axis transcript ID 90 degrees
      axis.text.y = element_text(size = 20, face = "bold.italic", angle = 90, hjust = 0.5, vjust = 0.5, color="black"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color="grey90", linetype="solid", linewidth = 0.5)
    )
}

# ==============================================================================
# 3. LOAD DATA AND PROCESS EXPRESSION
# ==============================================================================
cat("Loading and processing data...\n")

matrix_data <- read.table(opt$matrix, header = TRUE, check.names = FALSE)
colnames(matrix_data)[1] <- "gene_id"

matrix_data$gene_id <- case_when(
  matrix_data$gene_id == "NonamEVm000093t1" ~ "Ancestral CENP-E Paralog",
  matrix_data$gene_id == "NonamEVm002966t3" ~ "Truncated Neo-Y Transcript",
  TRUE ~ matrix_data$gene_id
)

metadata <- read.table(opt$meta, header = TRUE)

matrix_long <- matrix_data %>% pivot_longer(cols = -gene_id, names_to = "sample", values_to = "TMM")
df_combined <- inner_join(matrix_long, metadata, by = "sample")

df_summary <- df_combined %>%
  mutate(
    sex = case_when(
      sex %in% c("M", "Male", "male") ~ "Male",
      sex %in% c("F", "Female", "female") ~ "Female",
      TRUE ~ as.character(sex)
    ),
    tissue = case_when(
      tissue == "BR" ~ "Brain", tissue == "ES" ~ "Eyestalk", tissue == "GN" ~ "Gonad",
      tissue == "H"  ~ "Heart", tissue == "HP" ~ "Hepatopancreas", tissue == "MS" ~ "Muscle",
      tissue == "ST" ~ "Stomach", tissue == "WL" ~ "Walking Leg",
      TRUE ~ as.character(tissue)
    )
  ) %>%
  group_by(gene_id, tissue, sex) %>%
  summarise(
    mean_tmm = mean(TMM, na.rm = TRUE),
    se_tmm = ifelse(n() > 1, sd(TMM, na.rm = TRUE) / sqrt(n()), 0),
    .groups = "drop"
  )

# ==============================================================================
# 4. GENERATE ALL 4 PANELS
# ==============================================================================
cat("Generating sub-plots...\n")

# Top Row: Expression
plot_expr_ancestral <- generate_tmm_plot("Ancestral CENP-E Paralog", df_summary) + theme(axis.title.x = element_blank())
plot_expr_neoy <- generate_tmm_plot("Truncated Neo-Y Transcript", df_summary) + theme(axis.title.x = element_blank())

# Bottom Row: Protein Length 
length_ancestral <- 3020 
length_neoy      <- 780  

plot_len_ancestral <- generate_length_plot("NonamEVm000093t1", length_ancestral, len_colors["Ancestral"])
plot_len_neoy      <- generate_length_plot("NonamEVm002966t3", length_neoy, len_colors["NeoY"])

# ==============================================================================
# 5. COMBINE WITH PATCHWORK AND EXPORT
# ==============================================================================
cat("Stitching the final layout together...\n")

# Layout: Expression plots side-by-side, Length plots side-by-side underneath
top_row <- (plot_expr_ancestral | plot_expr_neoy)
bottom_row <- (plot_len_ancestral | plot_len_neoy)

# Combine and configure relative heights (Expression gets 80% of vertical space, Length gets 20%)
final_combined_figure <- (top_row / bottom_row) + 
  plot_layout(heights = c(4, 1), guides = "collect") & 
  theme(legend.position = "top")

# Define output filenames
png_out <- paste0(opt$out, ".png")
pdf_out <- paste0(opt$out, ".pdf")

# Save as a wide, high-resolution layout (24x14 inches)
ggsave(png_out, final_combined_figure, width=24, height=14, dpi=300, bg="white")
ggsave(pdf_out, final_combined_figure, width=24, height=14, device = cairo_pdf)

cat(paste("Success: 4-Panel combined figure saved as", png_out, "\n"))