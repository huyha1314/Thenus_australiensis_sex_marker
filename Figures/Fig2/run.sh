pixi run Rscript Fig2_A.R --input master_annotation_summary.tsv --out Fig2_A_Final
pixi run Rscript Fig2_B.R --input master_annotation_summary.tsv --out Fi2_B_Final
pixi run Rscript Fig2_combine.R --input master_annotation_summary.tsv --out Fig2_Combine_Final
