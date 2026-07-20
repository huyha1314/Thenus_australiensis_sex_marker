#!/usr/bin/env Rscript

# ==============================================================================
# Fig1S.R: SSR Visualization Tool for MISA Outputs
# ==============================================================================
# Description: This script parses MISA output files to generate premium, 
#              publication-quality genomics visualizations using ggplot2.
#              Follows the standardized project workflow structure.
# ==============================================================================

library(optparse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(scales)
library(gridExtra)
library(grid)

# Define command-line arguments using optparse
option_list <- list(
  make_option(c("-s", "--stats"), type = "character", default = "8.submit_scaffolds.fsa.statistics",
              help = "Path to MISA statistics file [default %default]", metavar = "character"),
  make_option(c("-m", "--misa"), type = "character", default = "8.submit_scaffolds.fsa.misa",
              help = "Path to raw MISA coordinate file [default %default]", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = "Fig1S_Final",
              help = "Output prefix for generated figures [default %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

stats_file <- opt$stats
misa_file <- opt$misa
out_prefix <- opt$out

cat("Configuration:\n")
cat("  Statistics file path:  ", stats_file, "\n")
cat("  MISA raw file path:    ", misa_file, "\n")
cat("  Output prefix:         ", out_prefix, "\n\n")

# Set elegant publication theme options globally
theme_set(
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 15)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "#555555", margin = margin(b = 20)),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10, color = "#333333"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#f0f0f0", linewidth = 0.5),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20)
  )
)

# Custom color palette (sleek and premium HSL-derived colors)
palette_ssr <- c(
  "Mono"     = "#4A90E2", # Deep Sky Blue
  "Di"       = "#50E3C2", # Mint Green
  "Tri"      = "#F5A623", # Warm Amber
  "Tetra"    = "#D0021B", # Vibrant Red
  "Penta"    = "#BD10E0", # Royal Purple
  "Hexa"     = "#9B9B9B", # Slate Gray
  "Compound" = "#7ED321"  # Lime Green
)

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------------------------
get_reverse_complement <- function(seq) {
  comp <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C", "N" = "N")
  seq_vec <- strsplit(toupper(seq), "")[[1]]
  comp_seq <- paste0(comp[seq_vec], collapse = "")
  comp_seq_rev <- paste(rev(strsplit(comp_seq, "")[[1]]), collapse = "")
  return(comp_seq_rev)
}

get_circular_permutations <- function(seq) {
  n <- nchar(seq)
  perms <- character(n)
  for (i in 1:n) {
    perms[i] <- paste0(substr(seq, i, n), substr(seq, 1, i - 1))
  }
  return(perms)
}

get_canonical_motif_class <- function(motif) {
  motif <- toupper(motif)
  perms_fwd <- get_circular_permutations(motif)
  perms_rev <- get_circular_permutations(get_reverse_complement(motif))
  
  all_candidates <- unique(c(perms_fwd, perms_rev))
  canonical_fwd <- min(perms_fwd)
  canonical_rev <- min(perms_rev)
  
  canonical_all <- min(all_candidates)
  partner <- if (canonical_all == canonical_fwd) canonical_rev else canonical_fwd
  partner <- min(get_circular_permutations(partner)) # Normalize partner
  
  if (canonical_all == partner) {
    return(paste0(canonical_all, "/", canonical_all))
  } else {
    return(paste0(canonical_all, "/", partner))
  }
}

# ------------------------------------------------------------------------------
# MISA FILE PARSERS
# ------------------------------------------------------------------------------
parse_misa_file <- function(filepath) {
  cat("Parsing raw MISA coordinate file:", filepath, "...\n")
  
  data <- read_tsv(filepath, col_types = cols(
    `ID` = col_character(),
    `SSR nr.` = col_double(),
    `SSR type` = col_character(),
    `SSR` = col_character(),
    `size` = col_double(),
    `start` = col_double(),
    `end` = col_double()
  ))
  
  colnames(data) <- c("scaffold", "ssr_num", "ssr_type", "ssr_string", "size", "start", "end")
  cat("  Loaded", nrow(data), "identified SSR records.\n")
  
  data <- data %>%
    mutate(
      class = case_when(
        ssr_type == "p1" ~ "Mono",
        ssr_type == "p2" ~ "Di",
        ssr_type == "p3" ~ "Tri",
        ssr_type == "p4" ~ "Tetra",
        ssr_type == "p5" ~ "Penta",
        ssr_type == "p6" ~ "Hexa",
        grepl("c", ssr_type) ~ "Compound",
        TRUE ~ "Compound"
      )
    )
  return(data)
}

# ------------------------------------------------------------------------------
# PLOTTING FUNCTIONS
# ------------------------------------------------------------------------------

# Plot 1: Pie Chart of SSR Class Distribution (Counts & Percentages)
plot_ssr_classes <- function(misa_data) {
  class_counts <- misa_data %>%
    group_by(class) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    mutate(
      Percentage = (Count / sum(Count)) * 100,
      class = factor(class, levels = c("Mono", "Di", "Tri", "Tetra", "Penta", "Hexa", "Compound"))
    ) %>%
    arrange(desc(class)) %>%
    mutate(
      label_pos = cumsum(Count) - (Count / 2)
    )
  
  total_ssrs <- sum(class_counts$Count)
  
  ggplot(class_counts, aes(x = "", y = Count, fill = class)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.6) +
    coord_polar("y", start = 0) +
    geom_text(aes(x = 1.2, y = label_pos, 
                  label = ifelse(Percentage > 2.0, 
                                 paste0(class, "\n", round(Percentage, 1), "%"), 
                                 "")), 
              color = "#111111", fontface = "bold", size = 3.5, lineheight = 0.9) +
    scale_fill_manual(values = palette_ssr, name = "SSR Class") +
    labs(
      title = "Distribution of SSR Classes",
      subtitle = paste("Total Identified SSRs =", comma(total_ssrs)),
      x = NULL,
      y = NULL
    ) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(t = 15, b = 5)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "#555555", margin = margin(b = 15)),
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      legend.position = "right",
      plot.margin = margin(15, 15, 15, 15)
    )
}

# ------------------------------------------------------------------------------
# MAIN PIPELINE EXECUTION
# ------------------------------------------------------------------------------
if (!file.exists(misa_file)) {
  stop(paste("Error: MISA raw coordinate file not found at:", misa_file))
}

misa_data <- parse_misa_file(misa_file)
p <- plot_ssr_classes(misa_data)

# Export sections
cat("\nSaving plots...\n")

pdf_output <- paste0(out_prefix, "_classes.pdf")
png_output <- paste0(out_prefix, "_classes.png")

# Save as vector PDF
ggsave(pdf_output, 
       plot = p, 
       width = 8, height = 7, 
       device = cairo_pdf)

# Save as high-resolution PNG
ggsave(png_output, 
       plot = p, 
       width = 8, height = 7, 
       dpi = 300, 
       bg = "white")

cat(paste("Success: Figures saved as", pdf_output, "and", png_output, "\n"))
