# Publication Figures & Reproducibility

This directory contains the scripts, data, and outputs required to generate the publication-ready figures for the *Thenus australiensis* sex marker discovery manuscript. 

We prioritize reproducibility and have built an automated, command-line-driven plotting pipeline using [Pixi](https://pixi.sh/) for cross-platform environment management and R scripts configured with `optparse`.

## Directory Structure
- **`Fig*`**: Main manuscript figures. Folders contain the R plotting scripts (e.g., `FigX.R`), necessary underlying data, and the execution script (`run.sh`).
- **`Fig*S`**: Supplementary figures following the same structure.
- **`Graphical_abstract`**: Final graphical abstract assets.

## 🛠 Prerequisites & Setup

To ensure you are using the exact same R version and package dependencies used in our study, you must first install the `pixi` package manager. 

**1. Install Pixi (if you haven't already)**
```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

**2. Initialize the plotting environment**
Navigate to this `Figures` directory and install the environment. This will automatically download and install R, `ggplot2`, `ggtree`, `patchwork`, and all other required bioconductor packages pinned in the `pixi.lock` file.
```bash
# Navigate to the Figures directory
cd /home/huyha/Desktop/lob/Thenus_australiensis_sex_marker/Figures

# Install the exact reproducible environment
pixi install
```

## 📊 How to Reproduce Figures

Each figure that has a reproducible R plotting script has been given its own directory. To reproduce a figure, simply navigate into its directory and execute the `run.sh` script. The `run.sh` script automatically wraps the R scripts in the `pixi run` command.

### Example: Reproducing Figure 3
```bash
# 1. Enter the specific figure's directory
cd Fig3

# 2. Execute the reproducible run script
bash run.sh
```

*This will execute the underlying R scripts and generate the output `Fig3_Final.pdf` and `Fig3_Final.png` directly inside the `Fig3` directory.*

### Running all figures
If you want to generate all figures that have execution scripts at once, you can run the following command from this directory:
```bash
for d in Fig*/; do if [ -f "$d/run.sh" ]; then echo "Running $d..."; (cd "$d" && bash run.sh); fi; done
```

## 📝 Modifying the Plots

If you wish to modify aesthetics (colors, dimensions, themes), you can directly edit the R scripts inside each folder. They are modularized using `optparse`, meaning you can also pass different inputs if you want to reuse our plotting architecture for your own datasets.

For example, checking the arguments for a script:
```bash
pixi run Rscript Fig3/Fig3.R --help
```
