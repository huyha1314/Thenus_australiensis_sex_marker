#!/bin/bash
source config.sh

echo "Starting 9. Genome assessment"
pixi run busco -c 10 \
  -m protein \
  --offline \
  --download_path busco_downloads \
  -i Thenus_australiensis_protein.faa \
  -l arthropoda_odb12 \
  -o busco_augustus
