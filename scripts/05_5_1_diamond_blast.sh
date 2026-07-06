#!/bin/bash
source config.sh

echo "Starting 5.1 Diamond Blast (v0.88.2) the whole NCBI NR database (Release 2021_9_29)"
pixi run diamond blastp \
  --query  anno/lobster_pasa_updated.proteins.no_stop.fasta \
  --db  db/nr.dmnd \
  --out nr.tsv \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
  --evalue 1e-5 \
  --max-target-seqs 1 \
  --threads 32
