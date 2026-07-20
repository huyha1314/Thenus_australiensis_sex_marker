#!/bin/bash
pixi run Rscript Fig3.R \
  --tree merged.aln.treefile \
  --ann Crustacea.tsv \
  --out Fig3_Final
