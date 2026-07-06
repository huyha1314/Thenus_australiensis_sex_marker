#!/bin/bash
source config.sh

echo "Starting 5.4 Protein annotation by eggNOG-Mappper"
pixi run emapper.py \
  -i  anno/gff/lobster_pasa_updated.proteins.fasta \
  --output  anno/agat_final/emmpper.out \
  -m diamond \
  --cpu 32 \
  --data_dir  db/eggnog_db
