#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)
library(patchwork)
library(zoo)
library(optparse)

# ==============================================================================
# 0. COMMAND LINE ARGUMENTS
# ==============================================================================
option_list = list(
  make_option(c("--wgsinput"), type="character", default="7345028_wgs_depth.txt", 
              help="Path to the WGS depth text file"),
  make_option(c("--dartinput"), type="character", default="7345028_dartseq_depth.txt", 
              help="Path to the DArT-seq depth text file"),
  make_option(c("--blastinput"), type="character", default="blast_result.txt", 
              help="Path to the BLAST result file"),
  make_option(c("--out"), type="character", default="Fig4_Updated_Final", 
              help="Prefix for the output files")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# ==============================================================================
# 1. PARSE BLAST RESULTS DYNAMICALLY
# ==============================================================================
blast_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                "qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast_data <- read.table(opt$blastinput, header=FALSE, col.names=blast_cols)

# Extract Neo-Y (Scc_2453) exons
neoy_blast <- blast_data %>% 
  filter(sseqid == "Scc_2453") %>%
  mutate(start = pmin(sstart, send), end = pmax(sstart, send), type = "Neo-Y Exon") %>%
  select(start, end, type)

# Extract Host (Scc_1163) exons
host_blast <- blast_data %>% 
  filter(sseqid == "Scc_1163") %>%
  mutate(start = pmin(sstart, send), end = pmax(sstart, send), type = "Host Exon") %>%
  select(start, end, type)

# Find max extents for the intron line
neoy_min <- min(neoy_blast$start, na.rm=TRUE)
neoy_max <- max(neoy_blast$end, na.rm=TRUE)
host_min <- min(host_blast$start, na.rm=TRUE)
host_max <- max(host_blast$end, na.rm=TRUE)

# ==============================================================================
# 2. BUILD PANEL A: THE OVERVIEW ARCHITECTURE
# ==============================================================================
host_navy    <- "#2A4B7C"
neo_y_green  <- "#009E73"
transposon_purple <- "#8E44AD"

transposons <- data.frame(
  x = c(58813, (58813+61548)/2, 61548, (58813+61548)/2, 104720, (104720+106467)/2, 106467, (104720+106467)/2),
  y = c(0.7, 0.85, 0.7, 0.55, 0.7, 0.85, 0.7, 0.55), 
  group = c(rep("Tc1",4), rep("LINE",4))
)

p_arch <- ggplot() +
  annotate("rect", xmin = 58813, xmax = 61548, ymin = 0.1, ymax = 0.9, fill = transposon_purple, alpha = 0.1) +
  annotate("text", x = (58813+61548)/2, y = 0.95, label = "Detailed Zoom Below", fontface = "italic", size = 6, color = transposon_purple) +
  
  # Neo-Y track
  geom_segment(aes(x = 0, xend = 122615, y = 0.7, yend = 0.7, color = "Intron"), linewidth = 1) +
  annotate("rect", xmin = neoy_min, xmax = neoy_max, ymin = 0.65, ymax = 0.75, fill = neo_y_green, alpha = 0.2) +
  geom_rect(data = neoy_blast, aes(xmin = start, xmax = end, ymin = 0.6, ymax = 0.8, fill = "Neo-Y Exon")) +
  geom_polygon(data = transposons, aes(x = x, y = y, group = group, fill = "Transposon")) +
  
  # Host track
  geom_segment(aes(x = 0, xend = host_max, y = 0.3, yend = 0.3, color = "Intron"), linewidth = 1) +
  annotate("rect", xmin = host_min, xmax = host_max, ymin = 0.25, ymax = 0.35, fill = host_navy, alpha = 0.2) +
  geom_rect(data = host_blast, aes(xmin = start, xmax = end, ymin = 0.2, ymax = 0.4, fill = "Host Exon")) +
  
  # Labels
  annotate("text", x = 2000, y = 0.83, label = "Truncated Neo-Y CENP-E (Expressed Testis mRNA - NonamEVm002966t3)", hjust = 0, fontface = "bold", color = neo_y_green, size = 6) +
  annotate("text", x = 2000, y = 0.43, label = "Ancestral Host CENP-E (Genomic DNA - NonamEVm000093t1)", hjust = 0, fontface = "bold", color = host_navy, size = 6) +
  
  annotate("text", x = (58813+61548)/2, y = 0.15, label = "Region A (Tc1) k141_7345028", color = transposon_purple, fontface = "bold", size = 6.25) +
  annotate("text", x = (104720+106467)/2, y = 0.15, label = "Region B (LINE/Jerky) k141_6893598", color = transposon_purple, fontface = "bold", size = 6.25) +
  
  coord_cartesian(xlim = c(0, 122615), ylim = c(0, 1)) + 
  labs(title = "A. Overview of Neo-Y Genomic Architecture (Scaffold_2453) vs. Host (Scaffold_1163)", x = "Genomic Position (bp)", y = "") +
  scale_fill_manual(name = "Feature Type:", values = c("Host Exon" = host_navy, "Neo-Y Exon" = neo_y_green, "Transposon" = transposon_purple)) +
  scale_color_manual(name = "", values = c("Intron" = "black")) +
  theme_classic() +
  theme(plot.title = element_text(size = 20, face = "bold"), axis.title.x = element_text(size = 19, face = "bold"), axis.text.x = element_text(size = 17), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom", legend.direction = "horizontal", legend.title = element_text(face = "bold", size = 21), legend.text = element_text(size = 20), plot.margin = margin(b = 12, unit = "pt")) 

# ==============================================================================
# 3. BUILD PANEL B: FACETED ZOOM PLOTS 
# ==============================================================================
wgs_data  <- read.table(opt$wgsinput, header = TRUE)
dart_data <- read.table(opt$dartinput, header = TRUE)

wgs_summary <- wgs_data %>%
  mutate(Mean_Female = rowMeans(select(., starts_with("F")), na.rm = TRUE),
         Mean_Male   = rowMeans(select(., starts_with("M")), na.rm = TRUE),
         DataType    = "WGS") %>% select(POS, Mean_Female, Mean_Male, DataType)

dart_summary <- dart_data %>%
  mutate(Mean_Female = rowMeans(select(., starts_with("Female")), na.rm = TRUE),
         Mean_Male   = rowMeans(select(., starts_with("Male")), na.rm = TRUE),
         DataType    = "DArT-seq") %>% select(POS, Mean_Female, Mean_Male, DataType)

combined_data <- bind_rows(wgs_summary, dart_summary)
smooth_window <- 10 

plot_data <- combined_data %>%
  group_by(DataType) %>%
  mutate(Smooth_Female = rollmean(Mean_Female, k = smooth_window, fill = "extend"),
         Smooth_Male   = rollmean(Mean_Male, k = smooth_window, fill = "extend")) %>% ungroup() %>%
  mutate(DataType = factor(DataType, levels = c("WGS", "DArT-seq")))

p_zoom <- ggplot(plot_data, aes(x = POS)) +
  annotate("rect", xmin = 1372, xmax = 1487, ymin = 0, ymax = Inf, fill = "grey80", alpha = 0.5) +
  geom_vline(xintercept = c(1372, 1487), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text(data = data.frame(POS = 1600, Smooth_Male = 25, DataType = "WGS"),
            aes(y = Smooth_Male, label = "PCR Target: 171 bp (1423-1593 bp)"), hjust = 0, fontface = "bold", size = 6) +
  geom_line(aes(y = Smooth_Male, color = "Male Cohort (Mean)"), linewidth = 1.2, alpha = 0.85) +
  geom_line(aes(y = Smooth_Female, color = "Female Cohort (Mean)"), linewidth = 1.2, alpha = 0.85) +
  facet_grid(DataType ~ ., scales = "free_y") +
  scale_color_manual(values = c("Male Cohort (Mean)" = "#1f77b4", "Female Cohort (Mean)" = "#D55E00")) +
  labs(title = "B. Detailed Zoom and PCR Target Validation (Contig k141_7345028)", x = "Position on Contig k141_7345028 (bp)", y = "Smoothed Read Depth", color = NULL) +
  theme_classic() +
  theme(strip.background = element_rect(fill = "grey95", color = "black"), strip.text = element_text(size = 19, face = "bold"), plot.title = element_text(size = 22, face = "bold"), axis.title.x = element_text(size = 19, face = "bold"), axis.title.y = element_text(size = 19, face = "bold"), axis.text = element_text(size = 17), legend.position = "bottom", legend.text = element_text(size = 22), panel.spacing = unit(1, "lines"), panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

combined_masterpiece <- p_arch / p_zoom + plot_layout(heights = c(1, 1.3))

pdf_output <- paste0(opt$out, ".pdf")
png_output <- paste0(opt$out, ".png")
ggsave(pdf_output, plot = combined_masterpiece, width = 14, height = 12, device = cairo_pdf)
ggsave(png_output, plot = combined_masterpiece, width = 14, height = 12, dpi = 300, bg = "white")
cat(paste("Success: Figures saved as", pdf_output, "and", png_output, "\n"))
