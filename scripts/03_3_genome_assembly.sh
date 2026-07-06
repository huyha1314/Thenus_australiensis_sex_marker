#!/bin/bash
source config.sh

echo "Starting 3. Genome assembly"
pixi run bbnorm.sh \
  in=merge_M4_R2.fastq.gz \
  out=/norm.merge_M4_R2.fastq.gz \
  target=40 mindepth=2 threads=24 -Xmx80G \
  hist=coverage.hist \
  prefilter=t passes=2

# Define input files
READS_R1="norm.merge_M4_R1.fastq.gz"
READS_R2="norm.merge_M4_R2.fastq.gz"
OUTPUT_DIR="abyss_results_65_M"
K_MER=65
# Use the number of CPUs requested

echo "Creating output directory: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

echo "Starting ABySS assembly with k=$K_MER and $NUM_PROCESSORS threads..."
echo "Input reads: $READS_R1, $READS_R2"

# Run abyss-pe
# -C . is redundant if you cd into the directory, but doesn't hurt.
# Removed H=3 to let ABySS determine optimal hash functions for B=30G
# Added v=-v for verbose output including Bloom filter stats
 pixi run abyss-pe \
  name=Male_65\
  k=$K_MER \
  in="$READS_R1 $READS_R2" \
  B=50G \
  kc=3 \
  v=-v \
  j=43 np=2
pixi run nextflow run nf-core/denovotranscript \
-r 1.2.1 \
--input samplesheet.csv \
--outdir $nf \
-profile singularity \
--assemblers trinity \
--ss fr \
--busco_mode transcriptome \
--busco_lineage arthropoda_odb10 \
-c run.config \
pixi run cd-hit-est -i nf/results/evigene/okayset/all_assembled.okay.mrna -o transcripts.cdhit97.fa \
  -c 0.97 -aS 0.95 -G 0 -g 1 -T 64 -M 0

pixi run salmon index -t transcripts.cdhit97.fa -i salmon/idx
pixi run salmon quant -i salmon/idx -l A \
-1 nf/results/cat/pooled_reads_1.merged.fastq.gz -2 nf/results/cat/pooled_reads_2.merged.fastq.gz \
-p 64 -o salmon

pixi run perl $TRINITY_HOME/util/abundance_estimates_to_matrix.pl \
  --est_method salmon \
  --name_sample_by_basedir \
  --quant_files salmon/quant_files.txt \
  --out_prefix matrix \
  --gene_trans_map none

# Keep transcripts with TPM ≥ 1 in at least one sample
pixi run perl $TRINITY_HOME/util/filter_low_expr_transcripts.pl \
  --matrix matrix.isoform.TMM.EXPR.matrix --transcripts busco/transcripts.cdhit97.fa --min_expr_any 1 \
  > salmon/okay.TPM1.fa

# Use a *stricter* identity (0.995) but keep 95% coverage of the shorter
pixi run cd-hit-est -i salmon/okay.TPM1.fa -o salmon/2okay.minredund.fa \
  -c 0.97 -aS 0.95 -G 0 -g 1 -T 64 -M 0
pixi run hisat2 -p 48 -x rnaseq/assembly -1 nf/results/cat/pooled_reads_1.merged.fastq.gz \
 -2  RNA/nf/results/cat/pooled_reads_2.merged.fastq.gz -S rnaseq/rna_aln_sam 
pixi run samtools sort -@25 -m 2G -o rnaseq/sort_rna_aln.bam rnaseq/rna_aln_sam
pixi run samtools index rnaseq/sort_rna_aln.bam

# scaffold using Rascaf
pixi run rascaf -b rnaseq/sort_rna_aln.bam -f rnaseq/rn_rascaf_scaffolded.fa -o rnaseq/rascaf_scaffolded.fa 
pixi run rascaf rascaf-join -r rnaseq/rascaf_scaffolded.fa.out -o rnaseq/rascaf_scaffolded
pixi run blat -stepSize=11 -repMatch=2253 -minScore=20 -minIdentity=80 rnaseq/rascaf_scaffolded.fa \
 salmon/okay.TPM1.fa \
  ./transcripts2.psl

pixi run bash L_RNA_scaffolder/L_RNA_scaffolder.sh \
 -d .  \
  -j rnaseq3/rn_rascaf_scaffolded.fa \
  -i L_RNA_scaffolder/transcripts2.noheader.psl \
  -o scaffolded_genome3
pixi run run-ntedit polish --draft scaffolded_genome/L_RNA_scaffolder.fasta --reads   data/merge/merge_M4 -k 31 \
-t 64  -f --solid

pixi run run-ntedit polish --draft nt/rename.ntedit_k31_edited.fa --read   data/merge/merge_M4 -k 41 \
-t 64  -f

pixi run run-ntedit polish --draft nt2/ntedit_k41_edited.fa --reads   data/merge/merge_M4 -k 25 \
-t 64  -f 
pixi run seqkit seq -m 2000 input.fasta > filtered_output.fasta
