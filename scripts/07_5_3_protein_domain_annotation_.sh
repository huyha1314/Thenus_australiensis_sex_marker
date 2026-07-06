#!/bin/bash
source config.sh

echo "Starting 5.3 Protein domain annotation by interproscan (v5.52-86.0) and pfam_scan.pl"
pixi run interproscan.sh \
  -i transcriptome_annotate/lobster_pasa_updated.proteins.no_stop.fasta \
  -cpu 64 \
  -goterms \
  -iprlookup \
  -pa
