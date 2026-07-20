#!/bin/bash
pixi run Rscript Fig3.R \
  --tree merged_mock_boot.treefile \
  --ann Crustacea.tsv \
  --out Fig3_Final
