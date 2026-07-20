library(png)
library(grid)
library(ggplot2)
library(patchwork)

read_img <- function(path) {
    img <- readPNG(path)
    rasterGrob(img, interpolate = TRUE)
}

f4 <- read_img("Fig4_Final.png")
f6 <- read_img("Fig4C.png")

p_f4 <- ggplot() + annotation_custom(f4) + theme_void()
p_f6 <- ggplot() + annotation_custom(f6) + theme_void()

final_fig2 <- p_f4 / p_f6 + 
    plot_layout(heights = c(1.5, 1)) + 
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 20, face = "bold"))

ggsave("Fig4_NeoY_Architecture_Validation.png", final_fig2, width = 12, height = 16, dpi = 300, bg = "white")
ggsave("Fig4_NeoY_Architecture_Validation.pdf", final_fig2, width = 12, height = 16, bg = "white")
