#!/bin/bash
# run_all.sh - Execute all steps of the pipeline sequentially
set -e

echo "Starting Genomics Pipeline..."

pixi run step_01_qc
pixi run step_02_survey
pixi run step_03_assembly
pixi run step_03_4_repeat
pixi run step_04_gene_pred
pixi run step_05_1_blast
pixi run step_05_2_blast
pixi run step_05_3_interpro
pixi run step_05_4_eggnog
pixi run step_06_ncrna
pixi run step_07_sex_marker
pixi run step_09_assessment
pixi run step_10_phylo

echo "Pipeline finished successfully!"
