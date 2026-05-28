#!/bin/bash
source config.sh

echo "Starting 4.Gene structure prediction"
# 1. Build HISAT2 index
seqkit seq -m 2000 scaffolded_genome3/assembly.fa > scaffolded_genome3/assembly.2kb.masked.fa

hisat2-build \
  scaffolded_genome3/assembly.2kb.masked.fa \
  lobster_index

# 2. Align RNA-seq reads
micromamba run -n hisat2 hisat2 -p 48 -x lobster_index \
  -1  RNA/nf/results/cat/pooled_reads_1.merged.fastq.gz \
  -2  RNA/nf/results/cat/pooled_reads_2.merged.fastq.gz | \
samtools sort -@10 -m 3G -o rna_alignments.bam
samtools index rna_alignments.bam

# 3. Run BRAKER with masked genome + RNA evidence
export GENEMARK_PATH= script_new/annotate/augustus/GeneMark-ETP/bin/gmes
export BAMTOOLS_PATH=/home/huyha/micromamba/envs/augustus/bin

 braker.pl \
  --genome= script_new/scaffolded_genome3/assembly.2kb.masked.fa  \
  --bam=rna_alignments.bam \
  --cores 84 \
  --species=slipper_lobster \
  --gff3 
singularity exec \
  -B /mnt/12T \
  pasapipeline_latest.sif \
  /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
  -c alignAssembly.config \
  -C -R \
  -g scaffolded_genome3/assembly.2kb.masked.fa \
  -t transcriptome_annotate/okay.TPM1.fa.clean \
  -T -u transcriptome_annotate/okay.TPM1.fa \
  --TDN tdn.accs \
  --ALIGNERS gmap,blat \
  --CPU 48
EVidenceModeler \
  --sample_id lobster4 \
  --genome assembly.2kb.masked.fa \
  --weights evm_weights.txt \
  --gene_predictions augustus.gff3 \
  --transcript_alignments sql.pasa.lite.pasa_assemblies.gff3 \
  --min_intron_length 20 \
  --segmentSize 1000000 \
  --overlapSize 10000 \
  --CPU 48

export PASAHOME=/home/huyha/micromamba/envs/pasa_env/opt/pasa-2.5.3

singularity exec \
  -B /mnt/12T \
  pasapipeline_latest.sif \
  /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
    -c alignAssembly.config \
    -A \
    -g scaffolded_genome3/assembly.2kb.masked.fa  \
    -t transcriptome_annotate/okay.TPM1.fa.clean \
    -L \
    --annots evidenModuler/lobster3.EVM.gff3 \
    --CPU 48

