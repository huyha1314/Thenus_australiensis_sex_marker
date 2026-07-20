#!/usr/bin/env bash
set -e
# Generate the annotation panel
pixi run Rscript Fig2_combine.R --input master_annotation_summary.tsv --out Fig2_Combine_Final
# Stitch with GenomeScope profiles
pixi run Rscript Fig2_stitch.R
