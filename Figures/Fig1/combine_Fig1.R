library(png)
library(grid)
library(ggplot2)
library(patchwork)
library(optparse)

# Parse command line options if any
option_list <- list(
  make_option(c("-o", "--out"), type="character", default="Fig1_Final", 
              help="Prefix for output files [default= %default]", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

read_img <- function(path) {
    img <- readPNG(path)
    rasterGrob(img, interpolate = TRUE)
}

# Read Panel A and Panel B
f1a <- read_img("1A.png")
f1b <- read_img("1B.png")

p_a <- ggplot() + annotation_custom(f1a) + theme_void()
p_b <- ggplot() + annotation_custom(f1b) + theme_void()

# Combine vertically (A on top, B on bottom)
final_fig <- (p_a / p_b) + 
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 24, face = "bold"))

# Save to PDF and PNG with a taller aspect ratio
png_out <- paste0(opt$out, ".png")
pdf_out <- paste0(opt$out, ".pdf")

ggsave(png_out, final_fig, width = 10, height = 14, dpi = 300, bg = "white")
ggsave(pdf_out, final_fig, width = 10, height = 14, bg = "white")

cat("Successfully combined Panel A and Panel B into", png_out, "and", pdf_out, "\n")
