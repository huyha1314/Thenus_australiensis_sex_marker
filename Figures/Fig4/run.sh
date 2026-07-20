#!/usr/bin/env bash
set -e
# Generate the genomic architecture plots
pixi run Rscript Fig4.R --wgsinput 7345028_wgs_depth.txt --dartinput 7345028_dartseq_depth.txt --out Fig4_Final
# Stitch with the gel image
pixi run Rscript Fig4_stitch.R
