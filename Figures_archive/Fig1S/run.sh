#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -e

# Run the Figure 1S script with standardized arguments
pixi run Rscript Fig1S.R \
  --stats "8.submit_scaffolds.fsa.statistics" \
  --misa "8.submit_scaffolds.fsa.misa" \
  --out "Fig1S_Final"
