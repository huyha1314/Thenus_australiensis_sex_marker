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
  make_option(c("--wgsinput"), type="character", default=NULL, 
              help="Path to the WGS depth text file", metavar="character"),
  make_option(c("--dartinput"), type="character", default=NULL, 
              help="Path to the DArT-seq depth text file", metavar="character"),
  make_option(c("--out"), type="character", default="Combined_Figure_Output", 
              help="Prefix for the output files [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$wgsinput) || is.null(opt$dartinput)){
  print_help(opt_parser)
  stop("Error: Both --wgsinput and --dartinput must be provided.", call.=FALSE)
}

wgs_input_file  <- opt$wgsinput
dart_input_file <- opt$dartinput
output_prefix   <- opt$out

# ==============================================================================
# 1. BUILD PANEL A: THE OVERVIEW ARCHITECTURE
# ==============================================================================
gene_start <- 4826
gene_end <- 122615
transcript_start <- 4841
transcript_end <- 86286

# Distinct color palette for Panel A
host_navy    <- "#2A4B7C"   # Navy Blue for Host Exons
neo_y_green  <- "#009E73"   # Green for mRNA
transposon_purple <- "#8E44AD" # Purple for Transposons

exons <- data.frame(start = c(4826, 10000, 20000, 35000, 45000, 65000, 75000, 90000, 115000, 122000),
                    end = c(5500, 10500, 21000, 36000, 46000, 66000, 76000, 91000, 116000, 122615), type = "Exon")
transcript_exons <- data.frame(start = c(4841, 10528, 13757, 17006, 23133, 35621, 40650, 42944, 53088, 59727, 67104, 75905, 84950, 85272, 85038, 85626, 85044),
                               end = c(4890, 10673, 13882, 17129, 23300, 35718, 40855, 43074, 53225, 59935, 67226, 75996, 85744, 85976, 85716, 86286, 85396), type = "RNA")
transposons <- data.frame(x = c(58813, (58813+61548)/2, 61548, (58813+61548)/2, 104720, (104720+106467)/2, 106467, (104720+106467)/2),
                          y = c(0.3, 0.45, 0.3, 0.15, 0.3, 0.45, 0.3, 0.15), group = c(rep("Tc1",4), rep("LINE",4)))

p_arch <- ggplot() +
  annotate("rect", xmin = 58813, xmax = 61548, ymin = 0.1, ymax = 0.9, fill = transposon_purple, alpha = 0.1) +
  annotate("text", x = (58813+61548)/2, y = 0.95, label = "Detailed Zoom Below", fontface = "italic", size = 6, color = transposon_purple) +
  
  geom_segment(aes(x = gene_start, xend = gene_end, y = 0.3, yend = 0.3, color = "Intron"), linewidth = 1) +
  geom_rect(data = exons, aes(xmin = start, xmax = end, ymin = 0.2, ymax = 0.4, fill = "Host Exon")) +
  geom_polygon(data = transposons, aes(x = x, y = y, group = group, fill = "Transposon")) +
  
  # FIXED: Mapped the Neo-Y backbone color to "Intron" so it inherits the black styling from the legend
  geom_segment(aes(x = transcript_start, xend = transcript_end, y = 0.7, yend = 0.7, color = "Intron"), linewidth = 0.8) +
  geom_rect(data = transcript_exons, aes(xmin = start, xmax = end, ymin = 0.6, ymax = 0.8, fill = "Neo-Y mRNA")) +
  
  annotate("text", x = gene_start, y = 0.45, label = "Ancestral Host CENP-E (Genomic DNA - NonamEVm000093t1)", hjust = 0, fontface = "bold", color = host_navy, size = 7) +
  annotate("text", x = (58813+61548)/2, y = 0.05, label = "Region A (Tc1) k141_7345028", color = transposon_purple, fontface = "bold", size = 6.25) +
  annotate("text", x = (104710+106467)/2, y = 0.05, label = "Region B (LINE/Jerky) k141_6893598", color = transposon_purple, fontface = "bold", size = 6.25) +
  annotate("text", x = 40000, y = 0.85, label = "Truncated Neo-Y CENP-E (Expressed Testis mRNA - NonamEVm002966t3)", hjust = 0.5, fontface = "bold", color = neo_y_green, size = 6) +
  
  coord_cartesian(xlim = c(0, 122615), ylim = c(0, 1)) + 
  labs(title = "A. Overview of Neo-Y Genomic Architecture (Scaffold_2453)", x = "Genomic Position on Scaffold_2453 (bp)", y = "") +
  
  scale_fill_manual(name = "Feature Type:", 
                    values = c("Host Exon" = host_navy, "Transposon" = transposon_purple, "Neo-Y mRNA" = neo_y_green)) +
  scale_color_manual(name = "", values = c("Intron" = "black")) +
  
  theme_classic() +
  theme(plot.title = element_text(size = 24, face = "bold"), 
        axis.title.x = element_text(size = 19, face = "bold"),
        axis.text.x = element_text(size = 17),
        axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom",          
        legend.direction = "horizontal",
        legend.title = element_text(face = "bold", size = 21),
        legend.text = element_text(size = 20),
        plot.margin = margin(b = 12, unit = "pt")) 

# ==============================================================================
# 2. BUILD PANEL B: FACETED ZOOM PLOTS 
# ==============================================================================
# Load data dynamically from the arguments
wgs_data  <- read.table(wgs_input_file, header = TRUE)
dart_data <- read.table(dart_input_file, header = TRUE)

# Calculate means and tag with DataType
wgs_summary <- wgs_data %>%
  mutate(Mean_Female = rowMeans(select(., starts_with("F")), na.rm = TRUE),
         Mean_Male   = rowMeans(select(., starts_with("M")), na.rm = TRUE),
         DataType    = "WGS") %>%
  select(POS, Mean_Female, Mean_Male, DataType)

dart_summary <- dart_data %>%
  mutate(Mean_Female = rowMeans(select(., starts_with("Female")), na.rm = TRUE),
         Mean_Male   = rowMeans(select(., starts_with("Male")), na.rm = TRUE),
         DataType    = "DArT-seq") %>%
  select(POS, Mean_Female, Mean_Male, DataType)

# Combine and smooth
combined_data <- bind_rows(wgs_summary, dart_summary)
smooth_window <- 10 

plot_data <- combined_data %>%
  group_by(DataType) %>%
  mutate(Smooth_Female = rollmean(Mean_Female, k = smooth_window, fill = "extend"),
         Smooth_Male   = rollmean(Mean_Male, k = smooth_window, fill = "extend")) %>%
  ungroup()

plot_data <- plot_data %>%
  mutate(DataType = factor(DataType, levels = c("WGS", "DArT-seq")))

# Create the Faceted Plot
p_zoom <- ggplot(plot_data, aes(x = POS)) +
  annotate("rect", xmin = 1372, xmax = 1487, ymin = 0, ymax = Inf, fill = "grey80", alpha = 0.5) +
  geom_vline(xintercept = c(1372, 1487), linetype = "dashed", color = "black", linewidth = 0.5) +
  
  geom_text(data = data.frame(POS = 1600, Smooth_Male = 25, DataType = "WGS"),
            aes(y = Smooth_Male, label = "PCR Target: 171 bp (1423-1593 bp)"), hjust = 0, fontface = "bold", size = 6) +
  
  geom_line(aes(y = Smooth_Male, color = "Male Cohort (Mean)"), linewidth = 1.2, alpha = 0.85) +
  geom_line(aes(y = Smooth_Female, color = "Female Cohort (Mean)"), linewidth = 1.2, alpha = 0.85) +
  
  facet_grid(DataType ~ ., scales = "free_y") +
  
  scale_color_manual(values = c("Male Cohort (Mean)" = "#1f77b4", "Female Cohort (Mean)" = "#D55E00")) +
  labs(title = "B. Detailed Zoom and PCR Target Validation (Contig k141_7345028)",
       x = "Position on Contig k141_7345028 (bp)", y = "Smoothed Read Depth", color = NULL) +
  
  theme_classic() +
  theme(strip.background = element_rect(fill = "grey95", color = "black"),
        strip.text = element_text(size = 19, face = "bold"),
        plot.title = element_text(size = 22, face = "bold"),
        axis.title.x = element_text(size = 19, face = "bold"),
        axis.title.y = element_text(size = 19, face = "bold"),
        axis.text = element_text(size = 17),
        legend.position = "bottom",
        legend.text = element_text(size = 22),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

# ==============================================================================
# 3. COMBINE EVERYTHING WITH PATCHWORK
# ==============================================================================

combined_masterpiece <- p_arch / p_zoom + plot_layout(heights = c(1, 1.3))

# --- EXPORT SECTION ---
pdf_output <- paste0(output_prefix, ".pdf")
png_output <- paste0(output_prefix, ".png")

# 1. Save as PDF
ggsave(pdf_output, 
       plot = combined_masterpiece, 
       width = 14, height = 12, 
       dpi = 300)

# 2. Save as PNG
ggsave(png_output, 
       plot = combined_masterpiece, 
       width = 14, height = 12, 
       dpi = 300, 
       bg = "white")

cat(paste("Success: Figures saved as", pdf_output, "and", png_output, "\n"))