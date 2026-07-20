#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(ggtree)
  library(ape)
  library(treeio)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(openxlsx)
  library(scales)
})

# ==============================================================================
# 0. COMMAND LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  make_option(c("-t", "--tree"), type="character", default=NULL, 
              help="Path to the tree file [e.g. merged.aln.treefile]", metavar="character"),
  make_option(c("-a", "--ann"), type="character", default=NULL, 
              help="Path to the annotation Crustacea.tsv file [e.g. Crustacea.tsv]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="Fig3_Final", 
              help="Prefix for the output files [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$tree) || is.null(opt$ann)){
  print_help(opt_parser)
  stop("Error: Both --tree and --ann must be provided.", call.=FALSE)
}

tree_input  <- opt$tree
ann_input   <- opt$ann
out_prefix  <- opt$out

# Verify files exist before proceeding
if (!file.exists(tree_input)) stop(paste("Error: Tree file not found at", tree_input))
if (!file.exists(ann_input)) stop(paste("Error: Annotation file not found at", ann_input))

# ==============================================================================
# 1. LOAD AND PREPROCESS DATA
# ==============================================================================
cat("Loading and processing phylogenetic tree and annotation metadata...\n")

# Load phylogenetic tree
tree <- read.tree(tree_input)

# Load annotation metadata (selecting labels and taxonomic info)
dat <- read.table(ann_input, header = TRUE, sep = "\t")[, c(2, 4:9)]
names(dat)[1] <- "label"

# Join tree metadata with phylogenetic tree
tree <- full_join(tree, dat, by = "label")
df <- na.omit(as_tibble(tree))

# ==============================================================================
# 2. BUILD THE PHYLOGENETIC TREE PLOT
# ==============================================================================
cat("Generating phylogenetic tree plot...\n")

# Set up the circular phylogenetic tree using ggtree
p <- ggtree(tree, layout = 'circular', linetype = 1, branch.length = 'none') + 
  geom_tree() +
  geom_point(aes(x + 0.5, y, shape = Order), size = 5) +
  geom_tiplab(aes(x + 1, y, 
                  label = Organism.Name, 
                  color = ifelse(label == "T.aus", "red", "black"),
                  fontface = ifelse(label == "T.aus", "bold", "plain")), 
              align = TRUE, size = 5) +
  geom_hilight(df, aes(node = node, fill = Class), type = "gradient", gradient.direction = 'rt', alpha = 0.8, to.bottom = TRUE) +
  scale_shape_manual(values = seq(0, 13)) +
  scale_color_identity() +
  xlim(0, 22) + # Adjusting the circular limits to give ample outer room for labels
  theme_tree(base_size = 15) +
  theme(
    legend.text = element_text(size = 15), 
    legend.title = element_text(size = 18, face = "bold"), 
    legend.key.size = unit(1.5, "cm"), 
    legend.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 8)), 
    fill  = guide_legend(override.aes = list(size = 8))
  )

# ==============================================================================
# 3. EXPORT SECTION
# ==============================================================================
pdf_output <- paste0(out_prefix, ".pdf")
png_output <- paste0(out_prefix, ".png")

cat("Saving plots...\n")

# 1. Save as PDF
ggsave(pdf_output, 
       plot = p, 
       width = 15, height = 15, 
       device = cairo_pdf)

# 2. Save as PNG (300 DPI with white background)
ggsave(png_output, 
       plot = p, 
       width = 15, height = 15, 
       dpi = 300, 
       bg = "white")

cat(paste("Success: Figures saved as", pdf_output, "and", png_output, "\n"))
