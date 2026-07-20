library(png)
library(grid)
library(ggplot2)
library(patchwork)

# Helper function to read and convert PNG to a ggplot-compatible grob
read_img <- function(path) {
    img <- readPNG(path)
    rasterGrob(img, interpolate = TRUE)
}

crop_white_borders <- function(path) {
    img <- readPNG(path)
    is_white <- apply(img, c(1, 2), function(pixel) {
        if (length(pixel) >= 4 && pixel[4] == 0) {
            return(TRUE)
        }
        all(pixel[1:3] > 0.99)
    })
    
    non_white_rows <- which(rowSums(!is_white) > 0)
    non_white_cols <- which(colSums(!is_white) > 0)
    
    if (length(non_white_rows) == 0 || length(non_white_cols) == 0) return(img)
    
    padding <- 10
    r_min <- max(1, min(non_white_rows) - padding)
    r_max <- min(nrow(img), max(non_white_rows) + padding)
    c_min <- max(1, min(non_white_cols) - padding)
    c_max <- min(ncol(img), max(non_white_cols) + padding)
    img[r_min:r_max, c_min:c_max, , drop = FALSE]
}

read_img_cropped <- function(path) {
    img <- crop_white_borders(path)
    rasterGrob(img, interpolate = TRUE)
}

cat("Running Update 3...\n")

f4 <- read_img("Fig4/Fig4_Final.png")
f6 <- read_img_cropped("Fig6/Fig6.png")

p_f4 <- wrap_elements(f4)
p_f6 <- wrap_elements(f6) + labs(tag = "C")

final_fig2 <- p_f4 / p_f6 + 
    plot_layout(heights = c(2.96, 1)) & 
    theme(plot.tag = element_text(size = 20, face = "bold"))

ggsave("New_Fig2_Genomic_Architecture_and_Gel.png", final_fig2, width = 12, height = 13.75, dpi = 300, bg = "white")
ggsave("New_Fig2_Genomic_Architecture_and_Gel.pdf", final_fig2, width = 12, height = 13.75, bg = "white")
ggsave("MER_Final_Lineup/Fig4_NeoY_Architecture_Validation.png", final_fig2, width = 12, height = 13.75, dpi = 300, bg = "white")
ggsave("MER_Final_Lineup/Fig4_NeoY_Architecture_Validation.pdf", final_fig2, width = 12, height = 13.75, bg = "white")

cat("Update 3 Complete!\n")
