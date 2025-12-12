library (VennDiagram)
library (dplyr)

# Read the clean list files
ids <- list(
  KEGG = scan("eggnog_ids.txt", what="", quiet=TRUE),
  NR = scan("nr_ids.txt", what="", quiet=TRUE),
  Swiss = scan("swiss_ids.txt", what="", quiet=TRUE),
  InterPro = scan("ipr_ids.txt", what="", quiet=TRUE)
)

# Plot
venn.diagram(
  x = ids,
  filename = "venn_diagram_output.png", # Saves directly to file
  col = c("blue4", "darkorchid4", "darkgoldenrod4", "antiquewhite4"),
  fill = c("turquoise", "pink", "gold", "bisque"),
  alpha = 0.45,
  cex = 1.25,
  cat.cex = 1.25,
  margin = 0.03
)