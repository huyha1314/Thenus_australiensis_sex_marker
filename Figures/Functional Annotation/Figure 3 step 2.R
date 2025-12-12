library (dplyr)
library (tidyr)
library (readr)
library (tidyverse)
library (ggplot2)
library(patchwork)

cog <- read.table ("/home/phuctran/Thenus_australiensis_sex_marker/COG_count.tsv")
colnames (cog) <- c('Category', 'Count')
cog_url <- "https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.fun.tab"
cog_key <- read.delim(cog_url,header = FALSE,sep = "\t",comment.char = "#",stringsAsFactors = FALSE)
colnames(cog_key) <- c("Category", "FunctionalGroup", "ColorHex", "Description")
merge <- merge (cog, cog_key, by ='Category')
merge <- merge %>%
  mutate(
    Group = case_when(
      Category %in% c("J", "A", "K", "L", "B") ~ "Information storage and processing",
      Category %in% c("D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O") ~ "Cellular processes and signaling",
      Category %in% c("C", "G", "E", "F", "H", "I", "P", "Q") ~ "Metabolism",
      Category %in% c("R", "S") ~ "Poorly characterized",
      TRUE ~ "Other"),
    CombinedLabel = paste0("[", Category, "] ", Description))

merge$Category <- factor(
  merge$Category, 
  levels = merge %>%
    arrange(Group, Count) %>% 
    pull(Category))

plot_data <- merge %>%
  arrange(Group, Count) %>%
  mutate(
    Category = factor(Category, levels = Category), 
    key_df = paste0("[", Category, "] ", Description))

p_plot <- ggplot(plot_data, aes(x = Category, y = Count, fill = Group)) +
  
  geom_bar(stat = "identity", color = "gray20", linewidth = 0.2) +
  geom_text(
    aes(label = Count), 
    vjust = -0.3, 
    color = "black",
    size = 2.5 ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  labs(
    title = "Functional Classification by COG",
    x = "COG Category Code", 
    y = "Number of Genes", 
    fill = "Major Functional Group") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.title.y = element_text(size = 10, margin = margin(r = 10)),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    panel.grid.major.y = element_blank())

key_data <- plot_data %>%
  arrange(Group, Category) %>%
  select(key_df) %>%
  pull()

num_rows <- ceiling(length(key_data) / 3)
col1 <- key_data[1:num_rows]
col2 <- key_data[(num_rows + 1):(2 * num_rows)]
col3 <- key_data[(2 * num_rows + 1):length(key_data)]


col2 <- c(col2, rep("", num_rows - length(col2)))
col3 <- c(col3, rep("", num_rows - length(col3)))

key_df_formatted <- data.frame(Col1 = col1, Col2 = col2, Col3 = col3)

y_coords <- 1:nrow(key_df_formatted)

p_key <- ggplot(key_df_formatted) +
  theme_void() +
  geom_text(aes(x = 0.02, y = y_coords, label = Col1), hjust = 0, size = 2.5) +
  geom_text(aes(x = 0.35, y = y_coords, label = Col2), hjust = 0, size = 2.5) +
  geom_text(aes(x = 0.68, y = y_coords, label = Col3), hjust = 0, size = 2.5) +
  xlim(0, 1) + 
  labs(title = "") +
  scale_y_reverse(limits = c(nrow(key_df_formatted) + 1, 0))

final_combined_figure <- p_plot / p_key + 
  plot_layout(heights = c(10, 3)) 

print(final_combined_figure)



