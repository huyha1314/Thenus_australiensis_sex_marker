#!/usr/bin/env Rscript

# Load required libraries
library(tidyverse)
library(optparse)

# ==============================================================================
# 0. COMMAND LINE ARGUMENTS
# ==============================================================================
option_list = list(
  make_option(c("--wgsinput"), type="character", default=NULL, 
              help="Path to the WGS depth text file", metavar="character"),
  make_option(c("--dartinput"), type="character", default=NULL, 
              help="Path to the DArT-seq depth text file", metavar="character"),
  make_option(c("--out"), type="character", default="Figure_Individual_Heatmap_Fixed", 
              help="Prefix for the output files [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$wgsinput) || is.null(opt$dartinput)){
  print_help(opt_parser)
  stop("Error: Both --wgsinput and --dartinput must be provided.", call.=FALSE)
}

# ====================================================================
# 1. Load and Format Data
# ====================================================================
wgs_data <- read.table(opt$wgsinput, header = TRUE)
wgs_long <- wgs_data %>%
  pivot_longer(cols = -c(CHROM, POS), names_to = "Sample", values_to = "Depth") %>%
  mutate(Technology = "WGS", Sex = ifelse(grepl("^F", Sample), "Female", "Male"))

dart_data <- read.table(opt$dartinput, header = TRUE)
dart_long <- dart_data %>%
  pivot_longer(cols = -c(CHROM, POS), names_to = "Sample", values_to = "Depth") %>%
  mutate(Technology = "DArT-seq", Sex = ifelse(grepl("Female", Sample), "Female", "Male"))

combined_data <- bind_rows(wgs_long, dart_long)

depth_cap <- 40
combined_data <- combined_data %>%
  mutate(Depth_Capped = ifelse(Depth > depth_cap, depth_cap, Depth),
         Technology = factor(Technology, levels = c("WGS", "DArT-seq")),
         Sex = factor(Sex, levels = c("Female", "Male")))

# ====================================================================
# 2. Generate the Fixed Heatmap Plot
# ====================================================================
marker_start <- 1372
marker_end   <- 1487

heatmap_plot <- ggplot(combined_data, aes(x = POS, y = Sample, fill = Depth_Capped)) +
  geom_raster() +
  scale_fill_gradient(low = "grey95", high = "#0072B2", 
                      name = "Read Depth",
                      limits = c(0, depth_cap),
                      na.value = "grey95") +
  geom_vline(xintercept = c(marker_start, marker_end), linetype = "dashed", color = "black", linewidth = 0.5) +
  
  # FIX: Using 'space = "free_y"' and 'scales = "free_y"' is correct, 
  # but we need to adjust the label formatting
  facet_grid(Technology + Sex ~ ., 
             scales = "free_y", 
             space = "free_y",
             labeller = labeller(Sex = c("Female" = "F", "Male" = "M"))) +
  
  labs(
    title = "Individual Sample Read Depth Mapping",
    subtitle = "Dashed lines indicate 115 bp male-specific region (1372-1487 bp)\nF = Female, M = Male",
    x = "Position on Contig k141_7345028 (bp)",
    y = "Individual Samples",
    fill = "Read Depth" 
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black"),
    # FIX: Increase the margin and adjust strip text size to prevent overlap
    strip.text.y = element_text(size = 16, face = "bold", angle = 270, margin = margin(l = 10, r = 10)),
    strip.placement = "outside",
    
    plot.title = element_text(size = 24, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 17, color = "black"),
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 17),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.spacing = unit(0.2, "lines") 
  )

# Export
ggsave(paste0(opt$out, ".png"), plot = heatmap_plot, width = 14, height = 10, dpi = 300, bg = "white")
ggsave(paste0(opt$out, ".pdf"), plot = heatmap_plot, width = 14, height = 10, dpi = 300, bg = "white")
