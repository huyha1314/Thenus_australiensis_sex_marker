library(png)
library(grid)
library(ggplot2)
library(patchwork)

read_img <- function(path) {
    img <- readPNG(path)
    rasterGrob(img, interpolate = TRUE)
}

f1_f <- read_img("F2_F.png")
f1_m <- read_img("F2_M.png")
f2 <- read_img("Fig2_Combine_Final_converted-1.png")

p1 <- ggplot() + annotation_custom(f1_f) + theme_void()
p2 <- ggplot() + annotation_custom(f1_m) + theme_void()
p3 <- ggplot() + annotation_custom(f2) + theme_void()

final_fig1 <- (p1 | p2) / p3 + 
    plot_layout(heights = c(1, 2.0)) + 
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 20, face = "bold"))

ggsave("Fig2_Genomic_Resource_Overview.png", final_fig1, width = 12, height = 18, dpi = 300, bg = "white")
ggsave("Fig2_Genomic_Resource_Overview.pdf", final_fig1, width = 12, height = 18, bg = "white")
