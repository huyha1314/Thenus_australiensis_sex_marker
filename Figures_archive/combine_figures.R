library(png)
library(grid)
library(ggplot2)
library(patchwork)

# Helper function to read and convert PNG to a ggplot-compatible grob
read_img <- function(path) {
    img <- readPNG(path)
    rasterGrob(img, interpolate = TRUE)
}

# Function to programmatically crop white space and transparent pixels from raster image arrays
crop_white_borders <- function(path) {
    img <- readPNG(path)
    # Identify pixels that are white (R, G, B > 0.99) or fully transparent (alpha == 0)
    is_white <- apply(img, c(1, 2), function(pixel) {
        if (length(pixel) >= 4 && pixel[4] == 0) {
            return(TRUE)
        }
        all(pixel[1:3] > 0.99)
    })
    
    non_white_rows <- which(rowSums(!is_white) > 0)
    non_white_cols <- which(colSums(!is_white) > 0)
    
    if (length(non_white_rows) == 0 || length(non_white_cols) == 0) {
        return(img)
    }
    
    # Crop the array with a small 10px padding to avoid clipping label borders
    padding <- 10
    r_min <- max(1, min(non_white_rows) - padding)
    r_max <- min(nrow(img), max(non_white_rows) + padding)
    c_min <- max(1, min(non_white_cols) - padding)
    c_max <- min(ncol(img), max(non_white_cols) + padding)
    
    img[r_min:r_max, c_min:c_max, , drop = FALSE]
}

# Helper to read and crop PNG to a ggplot-compatible grob
read_img_cropped <- function(path) {
    img <- crop_white_borders(path)
    rasterGrob(img, interpolate = TRUE)
}

# Update 2: Combine GenomeScope profiles (F1_F.png, F1_M.png) and functional annotation charts (Fig2_Combine_Final.png)
cat("Running Update 2...\n")

f1_f <- read_img("Fig1/F1_F.png")
f1_m <- read_img("Fig1/F1_M.png")
f2 <- read_img("Fig2/Fig2_Combine_Final.png")

p1 <- ggplot() + annotation_custom(f1_f) + theme_void()
p2 <- ggplot() + annotation_custom(f1_m) + theme_void()
p3 <- ggplot() + annotation_custom(f2) + theme_void()

# Combine using patchwork
# Horizontal top row (p1 | p2), then vertically stacked over p3
final_fig1 <- (p1 | p2) / p3 + 
    plot_layout(heights = c(1, 2.0)) + 
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 20, face = "bold"))

# Save to both paths to maintain backward compatibility and update target lineup files
ggsave("New_Fig1_Genomic_Resource_Overview.png", final_fig1, width = 12, height = 18, dpi = 300, bg = "white")
ggsave("New_Fig1_Genomic_Resource_Overview.pdf", final_fig1, width = 12, height = 18, bg = "white")
ggsave("MER_Final_Lineup/Fig2_Genomic_Resource_Overview.png", final_fig1, width = 12, height = 18, dpi = 300, bg = "white")
ggsave("MER_Final_Lineup/Fig2_Genomic_Resource_Overview.pdf", final_fig1, width = 12, height = 18, bg = "white")

cat("Update 2 Complete: Saved as Fig2_Genomic_Resource_Overview.png and .pdf in MER_Final_Lineup\n")


# Update 3: Merge Gel Image (Fig6.png) with Genomic Architecture (Fig4_Final.png)
cat("Running Update 3...\n")

f4 <- read_img("Fig4/Fig4_Final.png")
f6 <- read_img_cropped("Fig6/Fig6.png")

# Wrap elements to ensure full boundary preservation
p_f4 <- wrap_elements(f4)
p_f6 <- wrap_elements(f6) + labs(tag = "C")

# Combine using patchwork (stacked vertically)
# We calculated the mathematically optimal height ratio for zero vertical padding:
# H1 (Fig4) = W / 1.167 = 0.857 * W
# H2 (Fig6) = W / 3.46 = 0.289 * W
# Height ratio is 0.857 / 0.289 = 2.96 to 1.
# Total height = 12 * (0.857 + 0.289) = 13.75 inches.
final_fig2 <- p_f4 / p_f6 + 
    plot_layout(heights = c(2.96, 1)) & 
    theme(plot.tag = element_text(size = 20, face = "bold"))

# Save to both paths to maintain backward compatibility and update target lineup files
ggsave("New_Fig2_Genomic_Architecture_and_Gel.png", final_fig2, width = 12, height = 13.75, dpi = 300, bg = "white")
ggsave("New_Fig2_Genomic_Architecture_and_Gel.pdf", final_fig2, width = 12, height = 13.75, bg = "white")
ggsave("MER_Final_Lineup/Fig4_NeoY_Architecture_Validation.png", final_fig2, width = 12, height = 13.75, dpi = 300, bg = "white")
ggsave("MER_Final_Lineup/Fig4_NeoY_Architecture_Validation.pdf", final_fig2, width = 12, height = 13.75, bg = "white")

cat("Update 3 Complete: Saved as Fig4_NeoY_Architecture_Validation.png and .pdf in MER_Final_Lineup\n")

