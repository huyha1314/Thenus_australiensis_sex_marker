#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(VennDiagram))

# ==============================================================================
# 1. COMMAND LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type="character", default="master_annotation_summary.tsv", 
              help="Path to the master annotation TSV file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="Figure2_VennDiagram", 
              help="Prefix for the output file", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (!file.exists(opt$input)) {
  stop(paste("Error: Input file", opt$input, "not found."), call.=FALSE)
}

# ==============================================================================
# 2. READ DATA AND EXTRACT LISTS
# ==============================================================================
cat("Reading master annotation file...\n")
data <- read.delim(opt$input, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")

get_annotated_genes <- function(df, column_name) {
  df %>%
    filter(!!sym(column_name) != "-" & !!sym(column_name) != "") %>%
    pull(Gene_ID) %>%
    unique()
}

gene_lists <- list(
  KEGG      = get_annotated_genes(data, "KEGG_ko"),
  NR        = get_annotated_genes(data, "NR_Hit_ID"),
  SwissProt = get_annotated_genes(data, "SwissProt_Hit_ID"),
  InterPro  = get_annotated_genes(data, "InterPro_Acc")
)

# ==============================================================================
# 3. GENERATE THE VENN DIAGRAM 
# ==============================================================================
png_output <- paste0(opt$out, ".png")
pdf_output <- paste0(opt$out, ".pdf")

# Okabe-Ito Colorblind-Safe Palette (Sky Blue, Orange, Bluish-Green, Reddish-Purple)
cb_palette <- c("#56B4E9", "#E69F00", "#009E73", "#CC79A7")

pdf(NULL) # Prevent blank Rplots.pdf from generating

venn.plot <- venn.diagram(
  x = gene_lists,
  filename = NULL, 
  
  # Colors and Shapes 
  fill = cb_palette,
  alpha = 0.6,
  col = "black",         
  lwd = 1.5,             
  lty = "solid",
  
  # Number Styling (Inside the circles)
  cex = 1.5,             
  fontface = "bold",
  fontfamily = "sans",   
  
  # Category Label Styling (Outside the circles)
  cat.cex = 1.8,         
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  
  # FIX: Perfectly mapped to the 4 outer corners of the diagram
  # Order: KEGG (far-left), NR (top-right), SwissProt (top-left), InterPro (far-right)
  cat.pos = c(-150, 45, -45, 150), 
  
  # FIX: Pushing all 4 labels uniformly far away from the center
  cat.dist = c(0.26, 0.26, 0.26, 0.26), 
  
  margin = 0.15 
)

# ------------------------------------------------------------------------------
# 4. EXPORT TO PNG AND PDF
# ------------------------------------------------------------------------------

# Save as High-Res PNG
png(png_output, width = 3500, height = 3500, res = 300, bg = "white")
grid.draw(venn.plot)
dev.off()

# Save as High-Res PDF (Vector format)
pdf(pdf_output, width = 12, height = 12, bg = "white")
grid.draw(venn.plot)
dev.off()

cat(paste("Success: Publication-ready Venn diagrams saved as", png_output, "and", pdf_output, "\n"))