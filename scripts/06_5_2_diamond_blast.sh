#!/bin/bash
source config.sh

echo "Starting 5.2 Diamond Blast (v0.88.2) the SWISS-PROT database (Release 2022_01)"
diamond blastp \
  --query  anno/gff/lobster_pasa_updated.proteins.fasta  \
  --db  db/uniprot/uniprot_sprot.dmnd \
  --out swisprot.out \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
  --evalue 1e-5 \
  --max-target-seqs 1 \
  --threads 32
