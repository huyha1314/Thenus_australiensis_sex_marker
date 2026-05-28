#!/bin/bash
source config.sh

echo "Starting 6. Noncoding RNAs prediction"
seqkit split2 -p 20 -O chunks_dir scaffolded_genome3/rn_L_RNA_scaffolder.fasta

mkdir -p tRNAscan_chunks
mkdir -p $PWD/tmp_parallel
export TMPDIR=$PWD/tmp_parallel

parallel --tmpdir $PWD/tmp_parallel -j 20 '
  base=$(basename {} .fa);
  echo "Running tRNAscan-SE on $base";
  micromamba run -n rna_annot tRNAscan-SE --tmpdir rna/tmp  --thread 3 -m tRNAscan_chunks/${base}.stat \
              -o tRNAscan_chunks/${base}.tRNA \
              -f tRNAscan_chunks/${base}.structure \
              {}' ::: chunks_dir/*.fasta

cmsearch --cpu 24  --rfam --tblout genome.rfam.tbl \
  rna/db/Rfam.cm \
  scaffolded_genome3/rn_L_RNA_scaffolder.fasta > genome.fa.rfam
> 18S_ref.fasta
while read r; do
efetch -db nucleotide -id $r -format fasta >> 18S_ref.fasta
done < rna/snoscan/id.18s
> 28S_ref.fasta
while read n; do
efetch -db nucleotide -id $n -format fasta >> 28S_ref.fasta
done < rna/snoscan/id.28s
> ref_rRNA.fasta
cat 18S_ref.fasta 28S_ref.fasta > ref_rRNA.fasta
python3 make_snoscan_targets.py ref_rRNA.fasta targets.meth
seqkit split2 genome_wrapped.fasta -p 20 -O genome_parts

ls genome_parts/*.fasta | xargs -n 1 -P 20 -I {} sh -c 'snoscan -s ref_rRNA.fasta -m targets.meth "{}" > "{}.out"'

# 3. Merge results
cat genome_parts/*.out > snoscan_final.out


