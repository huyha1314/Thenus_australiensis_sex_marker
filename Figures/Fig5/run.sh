#!/usr/bin/env bash
set -e
pixi run Rscript Fig5.R --input align.txt --output Fig5_Transposon_Degradation_CENPE
