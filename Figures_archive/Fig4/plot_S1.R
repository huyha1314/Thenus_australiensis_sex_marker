#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)

data <- read.csv("/mnt/12T/lobster_project/script_new/test_sex_marker/2453/lobster_y_chromosome_candidates_separated.csv")

data_scaff <- data %>% filter(chrom == "scaffold_2453")

# We want to plot WGS depth across the scaffold to show the drop out
p1 <- ggplot(data_scaff, aes(x = start)) +
  # Highlight the transposon regions where drop-out occurs
  annotate("rect", xmin = 58813, xmax = 61548, ymin = -Inf, ymax = Inf, fill = "grey60", alpha = 0.3) +
  annotate("rect", xmin = 104720, xmax = 106467, ymin = -Inf, ymax = Inf, fill = "grey60", alpha = 0.3) +
  
  # Lines for Male and Female WGS depth
  geom_line(aes(y = wgs_m_median, color = "Male Cohort (Mean WGS Depth)"), linewidth = 1) +
  geom_line(aes(y = wgs_f_median, color = "Female Cohort (Mean WGS Depth)"), linewidth = 1, alpha = 0.7) +
  
  # Labels
  annotate("text", x = 60000, y = 3, label = "Region A (Tc1)\nMulti-mapping Drop-out", angle=90, fontface="bold") +
  annotate("text", x = 105500, y = 3, label = "Region B (LINE)\nMulti-mapping Drop-out", angle=90, fontface="bold") +
  
  labs(title = "Supplementary Figure S1: WGS Whole-Genome Mapping Profile (Scaffold_2453)",
       subtitle = "Standard MAPQ >= 30 filtering causes coverage drop-outs at highly repetitive transposon loci (Region A and B).",
       x = "Genomic Position (bp)", y = "Normalized Median Depth (WGS)", color = "") +
  
  theme_classic() +
  scale_color_manual(values = c("Male Cohort (Mean WGS Depth)" = "#1f77b4", "Female Cohort (Mean WGS Depth)" = "#D55E00")) +
  theme(legend.position="bottom", 
        legend.text = element_text(size=14),
        plot.title = element_text(size=18, face="bold"),
        plot.subtitle = element_text(size=14),
        axis.title = element_text(size=14, face="bold"),
        axis.text = element_text(size=12))

ggsave("Figure_S1.pdf", plot=p1, width=14, height=7, device=cairo_pdf)
ggsave("Figure_S1.png", plot=p1, width=14, height=7, dpi=300, bg="white")

cat("Successfully created Figure_S1.pdf and Figure_S1.png\n")
