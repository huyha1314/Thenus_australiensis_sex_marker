#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(stringr)) # Added for text wrapping
suppressPackageStartupMessages(library(optparse))

# ==============================================================================
# 1. COMMAND LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type="character", default="master_annotation_summary.tsv", 
              help="Path to the master annotation TSV file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="Figure5_COG_Classification", 
              help="Prefix for the output files", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (!file.exists(opt$input)) {
  stop(paste("Error: Input file", opt$input, "not found."), call.=FALSE)
}

# ==============================================================================
# 2. READ AND PROCESS COG DATA
# ==============================================================================
cat("Reading master annotation file and extracting COGs...\n")
data <- read.delim(opt$input, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")

# Extract the COG column, ignore hyphens and empties
cog_raw <- data %>%
  filter(COG_category != "-" & COG_category != "") %>%
  pull(COG_category)

# Split combined COGs (e.g., "GM" becomes "G" and "M")
cog_split <- unlist(strsplit(cog_raw, split = ""))

# Count the occurrences of each category
cog_counts <- data.frame(Category = cog_split) %>%
  group_by(Category) %>%
  summarise(Count = n(), .groups = 'drop')

# ==============================================================================
# 3. FETCH COG KEY AND MERGE
# ==============================================================================
cat("Fetching COG definitions from NCBI...\n")
cog_url <- "https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.fun.tab"
cog_key <- read.delim(cog_url, header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
colnames(cog_key) <- c("Category", "FunctionalGroup", "ColorHex", "Description")

# Merge data and assign major functional groups
merge_data <- merge(cog_counts, cog_key, by ='Category') %>%
  mutate(
    Group = case_when(
      Category %in% c("J", "A", "K", "L", "B") ~ "Information storage and processing",
      Category %in% c("D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O") ~ "Cellular processes and signaling",
      Category %in% c("C", "G", "E", "F", "H", "I", "P", "Q") ~ "Metabolism",
      Category %in% c("R", "S") ~ "Poorly characterized",
      TRUE ~ "Other"
    ),
    CombinedLabel = paste0("[", Category, "] ", Description)
  )

# Sort the data for the plot
plot_data <- merge_data %>%
  arrange(Group, Count) %>%
  mutate(
    Category = factor(Category, levels = Category), 
    key_df = paste0("[", Category, "] ", Description)
  )

# ==============================================================================
# 4. GENERATE THE BAR PLOT (PANEL A)
# ==============================================================================
p_plot <- ggplot(plot_data, aes(x = Category, y = Count, fill = Group)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.6) +
  
  # Text size on top of bars increased ~30% (from 5 to 6.5)
  geom_text(aes(label = Count), vjust = -0.5, color = "black", size = 6.5, fontface="bold") +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  scale_fill_brewer(palette = "Set2") + 
  
  labs(
    title = "COG Functional Classification",
    subtitle = "Thenus australiensis",
    x = "COG Category Code", 
    y = "Number of Genes", 
    fill = "Major Functional Group"
  ) +
  theme_classic() +
  
  # ALL TEXT SIZES INCREASED ~30%
  theme(
    plot.title = element_text(face = "bold", size = 31),
    plot.subtitle = element_text(size = 23),
    axis.title.x = element_text(size = 26, face = "bold", margin = margin(t=15)),
    axis.title.y = element_text(size = 26, face = "bold", margin = margin(r=15)),
    axis.text.x = element_text(size = 21, face = "bold", color = "black"),
    axis.text.y = element_text(size = 21, color = "black"),
    legend.title = element_text(size = 23, face = "bold"),
    legend.text = element_text(size = 21),
    legend.position = "top", 
    panel.grid.major.y = element_line(color="grey90", linetype="dashed")
  )

# ==============================================================================
# 5. GENERATE THE TEXT KEY (PANEL B) - OVERLAP FIX APPLIED
# ==============================================================================
# Extract the descriptions and WRAP the text at 45 characters
key_data <- plot_data %>%
  arrange(Group, Category) %>%
  pull(key_df) %>%
  str_wrap(width = 45) # This forces long lines to break into multiple lines

# Create 3 columns for the legend key
num_rows <- ceiling(length(key_data) / 3)
col1 <- key_data[1:num_rows]
col2 <- key_data[(num_rows + 1):(2 * num_rows)]
col3 <- key_data[(2 * num_rows + 1):length(key_data)]

# Pad with empty strings if columns are uneven
col2 <- c(col2, rep("", num_rows - length(col2)))
col3 <- c(col3, rep("", num_rows - length(col3)))

key_df_formatted <- data.frame(Col1 = col1, Col2 = col2, Col3 = col3)

# Increase vertical spacing so wrapped lines don't overlap vertically
y_coords <- seq(from = 1, by = 2.5, length.out = nrow(key_df_formatted))

# Legend text size increased ~30% (from 5 to 6.5)
p_key <- ggplot(key_df_formatted) +
  theme_void() +
  geom_text(aes(x = 0.01, y = y_coords, label = Col1), hjust = 0, vjust = 1, size = 6.5, lineheight = 0.8) +
  geom_text(aes(x = 0.34, y = y_coords, label = Col2), hjust = 0, vjust = 1, size = 6.5, lineheight = 0.8) +
  geom_text(aes(x = 0.67, y = y_coords, label = Col3), hjust = 0, vjust = 1, size = 6.5, lineheight = 0.8) +
  xlim(0, 1) + 
  scale_y_reverse(limits = c(max(y_coords) + 2, 0)) + # Expanded height to accommodate the wrapped text
  theme(plot.margin = margin(t = 20, b = 20))

# ==============================================================================
# 6. COMBINE AND EXPORT
# ==============================================================================
# Increased the bottom proportion to fit the newly wrapped text block
final_combined_figure <- p_plot / p_key + plot_layout(heights = c(10, 6)) 

pdf_output <- paste0(opt$out, ".pdf")
png_output <- paste0(opt$out, ".png")

# Increased canvas size (18x16) to gracefully hold the 30% larger fonts
ggsave(pdf_output, plot = final_combined_figure, width = 18, height = 16, dpi = 300)
ggsave(png_output, plot = final_combined_figure, width = 18, height = 16, dpi = 300, bg = "white")

cat(paste("Success: COG classification figures saved as", pdf_output, "and", png_output, "\n"))