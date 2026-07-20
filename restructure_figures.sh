#!/usr/bin/env bash
set -e

BASE_DIR="/home/huyha/Desktop/lob/Thenus_australiensis_sex_marker"
OLD_FIGS="${BASE_DIR}/Figures_archive"
NEW_FIGS="${BASE_DIR}/Figures"

echo "1. Archiving current Figures to Figures_archive..."
if [ ! -d "$OLD_FIGS" ]; then
    mv "$NEW_FIGS" "$OLD_FIGS"
else
    echo "Figures_archive already exists, proceeding..."
fi

echo "2. Re-creating new Figures directory..."
rm -rf "$NEW_FIGS"
mkdir -p "$NEW_FIGS"

echo "3. Copying environment & root files..."
cp "${OLD_FIGS}/pixi.toml" "${OLD_FIGS}/pixi.lock" "${OLD_FIGS}/.gitignore" "${OLD_FIGS}/.gitattributes" "${OLD_FIGS}/README.md" "$NEW_FIGS/"

echo "4. Building New Figure 1 (Bioinformatics Workflow)..."
mkdir -p "${NEW_FIGS}/Fig1"
cat << 'EOF' > "${NEW_FIGS}/Fig1/README.md"
# Figure 1: The Bioinformatic Workflow

**Status**: [PENDING MANUAL DESIGN]

This figure is the only one in the *Molecular Ecology Resources* (MER) lineup that cannot and should not be generated using R/Python code. Programmatic flowcharts are rigid and lack the required biological iconography.

**Instructions**:
1. Open a visual vector editor like **BioRender** or **Lucidchart**.
2. Create a 3-column flowchart: WGS, Transcriptomics (RNA-seq), and Genotyping (DArT-seq).
3. Draw arrows converging into a central step: **k-mer subtraction & structural mapping**.
4. End with a final "Validation" box highlighting the **171 bp marker**.
5. Save the exported files in this directory as `Fig1_Bioinformatics_Workflow.png` and `.pdf`.
EOF


echo "5. Building New Figure 2 (Genomic Resource Overview)..."
mkdir -p "${NEW_FIGS}/Fig2"
cp "${OLD_FIGS}/Fig1/"*.png "${NEW_FIGS}/Fig2/" || echo "No PNG in Fig1"
cp "${OLD_FIGS}/Fig2/Fig2_combine.R" "${OLD_FIGS}/Fig2/master_annotation_summary.tsv" "${NEW_FIGS}/Fig2/"

# Create the stitching script for Fig2
cat << 'EOF' > "${NEW_FIGS}/Fig2/Fig2_stitch.R"
library(png)
library(grid)
library(ggplot2)
library(patchwork)

read_img <- function(path) {
    img <- readPNG(path)
    rasterGrob(img, interpolate = TRUE)
}

# Adjusted to load 1A.png since F1_F.png and F1_M.png were deleted
f1 <- read_img("1A.png")
f2 <- read_img("Fig2_Combine_Final.png")

p1 <- ggplot() + annotation_custom(f1) + theme_void()
p3 <- ggplot() + annotation_custom(f2) + theme_void()

final_fig1 <- p1 / p3 + 
    plot_layout(heights = c(1, 2.0)) + 
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 20, face = "bold"))

ggsave("Fig2_Genomic_Resource_Overview.png", final_fig1, width = 12, height = 18, dpi = 300, bg = "white")
ggsave("Fig2_Genomic_Resource_Overview.pdf", final_fig1, width = 12, height = 18, bg = "white")
EOF

cat << 'EOF' > "${NEW_FIGS}/Fig2/run.sh"
#!/usr/bin/env bash
set -e
# Generate the annotation panel
pixi run Rscript Fig2_combine.R --input master_annotation_summary.tsv --out Fig2_Combine_Final
# Stitch with GenomeScope profiles
pixi run Rscript Fig2_stitch.R
EOF
chmod +x "${NEW_FIGS}/Fig2/run.sh"

echo "6. Building New Figure 3 (Phylogenetic Context)..."
cp -r "${OLD_FIGS}/Fig3" "${NEW_FIGS}/Fig3"

echo "7. Building New Figure 4 (Neo-Y Architecture)..."
mkdir -p "${NEW_FIGS}/Fig4"
cp "${OLD_FIGS}/Fig4/"* "${NEW_FIGS}/Fig4/"
cp "${OLD_FIGS}/Fig6/Fig6.png" "${NEW_FIGS}/Fig4/"

cat << 'EOF' > "${NEW_FIGS}/Fig4/Fig4_stitch.R"
library(png)
library(grid)
library(ggplot2)
library(patchwork)

read_img <- function(path) {
    img <- readPNG(path)
    rasterGrob(img, interpolate = TRUE)
}

f4 <- read_img("Fig4_Final.png")
f6 <- read_img("Fig6.png")

p_f4 <- ggplot() + annotation_custom(f4) + theme_void()
p_f6 <- ggplot() + annotation_custom(f6) + theme_void()

final_fig2 <- p_f4 / p_f6 + 
    plot_layout(heights = c(1.5, 1)) + 
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 20, face = "bold"))

ggsave("Fig4_NeoY_Architecture_Validation.png", final_fig2, width = 12, height = 16, dpi = 300, bg = "white")
ggsave("Fig4_NeoY_Architecture_Validation.pdf", final_fig2, width = 12, height = 16, bg = "white")
EOF

cat << 'EOF' > "${NEW_FIGS}/Fig4/run.sh"
#!/usr/bin/env bash
set -e
# Generate the genomic architecture plots
pixi run Rscript Fig4.R --wgsinput 7345028_wgs_depth.txt --dartinput 7345028_dartseq_depth.txt --out Fig4_Final
# Stitch with the gel image
pixi run Rscript Fig4_stitch.R
EOF
chmod +x "${NEW_FIGS}/Fig4/run.sh"

echo "8. Building New Figure 5 (CENP-E Degradation, was Fig3S)..."
mkdir -p "${NEW_FIGS}/Fig5"
cp "${OLD_FIGS}/Fig3S/"* "${NEW_FIGS}/Fig5/"
# Rename script to Fig5.R
if [ -f "${NEW_FIGS}/Fig5/Fig3S.R" ]; then
    mv "${NEW_FIGS}/Fig5/Fig3S.R" "${NEW_FIGS}/Fig5/Fig5.R"
fi
cat << 'EOF' > "${NEW_FIGS}/Fig5/run.sh"
#!/usr/bin/env bash
set -e
pixi run Rscript Fig5.R --input align.txt --output Fig5_Final
EOF
chmod +x "${NEW_FIGS}/Fig5/run.sh"

echo "9. Building New Figure 6 (Expression Divergence, was Fig5)..."
mkdir -p "${NEW_FIGS}/Fig6"
cp "${OLD_FIGS}/Fig5/"* "${NEW_FIGS}/Fig6/"
if [ -f "${NEW_FIGS}/Fig6/Fig5.R" ]; then
    mv "${NEW_FIGS}/Fig6/Fig5.R" "${NEW_FIGS}/Fig6/Fig6.R"
fi
cat << 'EOF' > "${NEW_FIGS}/Fig6/run.sh"
#!/usr/bin/env bash
set -e
pixi run Rscript Fig6.R
EOF
chmod +x "${NEW_FIGS}/Fig6/run.sh"

echo "10. Building New Figure 7 (Population Validation)..."
cp -r "${OLD_FIGS}/Fig7" "${NEW_FIGS}/Fig7"
# Fix the run.sh if it had issues, though we just copy it for now.

echo "11. Building Fig1S..."
cp -r "${OLD_FIGS}/Fig1S" "${NEW_FIGS}/Fig1S"

echo "Reorganization complete!"
