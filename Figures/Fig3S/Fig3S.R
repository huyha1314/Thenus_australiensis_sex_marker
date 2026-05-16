#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(ggplot2)
})

# 1. Define Command-Line Arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path to the input MSA file (e.g., align.txt)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="alignment_plot.png", 
              help="Path for the output plot image [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Input file must be supplied using --input.n", call.=FALSE)
}

# 2. Load the MSA file
cat("Loading alignment from:", opt$input, "\n")
aligned_seqs <- readDNAStringSet(opt$input)
names(aligned_seqs) <- sapply(strsplit(names(aligned_seqs), " "), `[`, 1)

# Extract sequences
char_long <- strsplit(as.character(aligned_seqs["NonamEVm000093t1"]), "")[[1]]
char_short <- strsplit(as.character(aligned_seqs["NonamEVm002966t3"]), "")[[1]]

# 3. Create the binary coverage data frame
align_df <- data.frame(
  Position = 1:length(char_long),
  Long_Cov = ifelse(char_long != "-", 1, 0),
  Short_Cov = ifelse(char_short != "-", 1, 0)
)
align_df$Overlap <- ifelse(align_df$Long_Cov == 1 & align_df$Short_Cov == 1, 1, 0)

# Function to find contiguous sequence blocks for clean rendering
get_blocks <- function(cov_vector, ymin, ymax) {
  rle_res <- rle(cov_vector)
  ends <- cumsum(rle_res$lengths)
  starts <- c(1, ends[-length(ends)] + 1)
  
  blocks <- data.frame(Start = starts, End = ends, Value = rle_res$values)
  blocks <- blocks[blocks$Value == 1, ]
  
  if(nrow(blocks) > 0) {
      blocks$ymin <- ymin
      blocks$ymax <- ymax
      return(blocks)
  } else {
      return(NULL)
  }
}

long_blocks <- get_blocks(align_df$Long_Cov, 2.0, 2.6)
short_blocks <- get_blocks(align_df$Short_Cov, 1.0, 1.6)
overlap_blocks <- get_blocks(align_df$Overlap, 0.0, 0.6)

# 4. Map eggNOG coordinates and make Domains "Gap-Aware"
long_non_gaps <- which(char_long != "-")
short_non_gaps <- which(char_short != "-")

# Specific coordinates from eggNOG output
long_msa_start <- long_non_gaps[333]
long_msa_end <- long_non_gaps[2216]
short_msa_start <- short_non_gaps[150]
short_msa_end <- short_non_gaps[1049]

# NEW: Create a true/false map of where the domains actually have sequence
align_df$Long_Domain <- ifelse(align_df$Position >= long_msa_start & align_df$Position <= long_msa_end & align_df$Long_Cov == 1, 1, 0)
align_df$Short_Domain <- ifelse(align_df$Position >= short_msa_start & align_df$Position <= short_msa_end & align_df$Short_Cov == 1, 1, 0)

# Extract only the solid blocks for the domains (ignores gaps)
long_domain_blocks <- get_blocks(align_df$Long_Domain, 2.1, 2.5)
short_domain_blocks <- get_blocks(align_df$Short_Domain, 1.1, 1.5)

# Label positioning
domain_labels_df <- data.frame(
  Transcript = c("Ancestral Motor Domain", "Truncated Motor Core"),
  X_pos = c(long_msa_start + (long_msa_end - long_msa_start)/2, 
            short_msa_start + (short_msa_end - short_msa_start)/2),
  Y_pos = c(2.65, 1.65)
)

# 5. Define Colorblind & Print-Safe Palette (Okabe-Ito based)
col_long    <- "#0072B2"  # Dark Blue
col_short   <- "#56B4E9"  # Light Blue
col_overlap <- "#000000"  # Black
col_domain  <- "#E69F00"  # Orange

# 6. Generate the Plot
cat("Generating plot...\n")
p <- ggplot() +
  # Sequence Blocks (Background)
  geom_rect(data = long_blocks, aes(xmin = Start, xmax = End, ymin = ymin, ymax = ymax),
            fill = col_long, color = NA) + 
  geom_rect(data = short_blocks, aes(xmin = Start, xmax = End, ymin = ymin, ymax = ymax),
            fill = col_short, color = NA) + 
  geom_rect(data = overlap_blocks, aes(xmin = Start, xmax = End, ymin = ymin, ymax = ymax),
            fill = col_overlap, color = NA) +
  
  # NEW: Gap-Aware Domains
  geom_rect(data = long_domain_blocks, aes(xmin = Start, xmax = End, ymin = ymin, ymax = ymax),
            fill = col_domain, color = col_overlap, linewidth = 0.2) +
  geom_rect(data = short_domain_blocks, aes(xmin = Start, xmax = End, ymin = ymin, ymax = ymax),
            fill = col_domain, color = col_overlap, linewidth = 0.2) +
            
  # Labels
  geom_text(data = domain_labels_df, aes(x = X_pos, y = Y_pos, label = Transcript),
            fontface = "bold", size = 4.5) +

  # Formatting & Labels
  scale_y_continuous(breaks = c(0.3, 1.3, 2.3), 
                     labels = c("Alignment Overlap", 
                                "Neo-Y CENP-E\n(NonamEVm002966t3)", 
                                "Ancestral CENP-E\n(NonamEVm000093t1)")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "CENP-E Transcript Alignment with eggNOG Motor Domains",
    subtitle = "Highlighting sequence divergence and structural truncation of the neo-Y paralog",
    x = "Alignment Position (bp)",
    y = ""
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )
  
# 7. Save the output in multiple formats
cat("Saving plots...\n")

# Extract the base name (removes .png or .pdf if you accidentally typed it)
base_name <- tools::file_path_sans_ext(opt$output)

# Save as PNG
png_file <- paste0(base_name, ".png")
cat(" -> Saving PNG to:", png_file, "\n")
ggsave(png_file, plot = p, width = 12, height = 6, dpi = 300, bg = "white")

# Save as PDF (PDFs are vector graphics, so DPI is not needed)
pdf_file <- paste0(base_name, ".pdf")
cat(" -> Saving PDF to:", pdf_file, "\n")
ggsave(pdf_file, plot = p, width = 12, height = 6, bg = "white")

cat("Done! Both formats generated successfully.\n")