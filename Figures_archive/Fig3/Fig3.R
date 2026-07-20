#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(ggtree)
  library(ape)
  library(treeio)
  library(dplyr)
  library(ggplot2)
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

tree <- read.tree(tree_input)
dat <- read.table(ann_input, header = TRUE, sep = "\t")[, c(2, 4:9)]
names(dat)[1] <- "label"

# Prepare label_expr before joining
dat$label_expr <- ifelse(
  dat$label == "T.aus",
  'italic("Thenus australiensis")~bold(" (This study)")',
  paste0('italic("', dat$Organism.Name, '")')
)

tree_joined <- full_join(tree, dat, by = "label")

# Mutate to add target and bootstrap columns
phylo_tree <- as.phylo(tree_joined)
n_tips <- Ntip(phylo_tree)
t_aus_idx <- which(phylo_tree$tip.label == "T.aus")
target_nodes <- nodepath(phylo_tree, from = n_tips + 1, to = t_aus_idx)

tree_joined <- tree_joined %>%
  mutate(
    is_target_branch = node %in% target_nodes,
    branch_color = ifelse(is_target_branch, "#D55E00", "grey40"),
    is_target_label = (label == "T.aus"),
    bootstrap = suppressWarnings(ifelse(node > n_tips & !is.na(label) & label != "", as.numeric(label), NA)),
    support_color = case_when(
      node <= n_tips ~ NA_character_,
      is.na(bootstrap) ~ NA_character_,
      bootstrap >= 95 ~ "black",
      bootstrap >= 80 ~ "grey60",
      TRUE ~ NA_character_
    )
  )

# Custom palette
custom_palette <- c(
  "Branchiopoda" = "#D55E00",
  "Hexanauplia"  = "#E69F00",
  "Malacostraca" = "#0072B2",
  "Ostracoda"    = "#56B4E9",
  "Thecostraca"  = "#CC79A7"
)

# Hilight data: filter out rows without HilightClass, and avoid dropping rows due to other NA columns
df_hilight <- as_tibble(tree_joined) %>%
  filter(!is.na(Class) & Class != "") %>%
  rename(HilightClass = Class)

# ==============================================================================
# 2. BUILD THE PHYLOGENETIC TREE PLOT
# ==============================================================================
cat("Generating phylogenetic tree plot...\n")

order_shapes <- c(
  "Amphipoda"         = 0,  # square
  "Anostraca"         = 1,  # circle
  "Balanomorpha"      = 2,  # triangle point up
  "Calanoida"         = 8,  # asterisk (distinct from Decapoda)
  "Decapoda"          = 15, # solid square (highly visible for the main focal group)
  "Diplostraca"       = 5,  # diamond
  "Euphausiacea"      = 6,  # triangle point down
  "Harpacticoida"     = 7,  # square with cross
  "Isopoda"           = 9,  # diamond with plus
  "Podocopida"        = 10, # circle with plus
  "Pollicipedomorpha" = 11, # star
  "Siphonostomatoida" = 12, # square with plus
  "Stomatopoda"       = 13  # circle with cross
)

p <- ggtree(tree_joined, layout = 'circular', linetype = 1, aes(color = branch_color, linewidth = is_target_branch)) + 
  scale_linewidth_manual(values = c("FALSE" = 0.7, "TRUE" = 2.0), guide = "none") +
  scale_color_identity() +
  geom_hilight(data = df_hilight, aes(node = node, fill = HilightClass), type = "gradient", gradient.direction = 'rt', alpha = 0.4, to.bottom = TRUE) +
  scale_fill_manual(name = "Taxonomic Class", values = custom_palette, na.translate = FALSE) +
  geom_nodepoint(aes(color = support_color), size = 2.0, na.rm = TRUE) +
  geom_point(data = function(x) subset(x, isTip), aes(x = 3.1, y = y, shape = Order), size = 4, stroke = 1.2, color = "black") +
  scale_shape_manual(values = order_shapes, na.translate = FALSE) +
  geom_tiplab(aes(x = 3.3, label = label_expr, color = ifelse(is_target_label, "red", "black")), 
              align = TRUE, size = 4, parse = TRUE, linetype = "dotted", linesize = 0.3) +
  geom_treescale(x = 0.8, y = 0, width = 0.2, fontsize = 5, linesize = 1.0, offset = 0.05) +
  xlim(0, 5.5) +
  theme_tree(base_size = 14) +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(1.0, "cm")
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
       width = 12, height = 12, 
       device = cairo_pdf)

# 2. Save as PNG (300 DPI with white background)
ggsave(png_output, 
       plot = p, 
       width = 12, height = 12, 
       dpi = 300, 
       bg = "white")

cat(paste("Success: Figures saved as", pdf_output, "and", png_output, "\n"))
