#!/bin/bash
source config.sh

echo "Starting Repeat prediction and annotation"
# run earlygrey build database 
pixi run earlGrey \
  -g /worker_data/huyha/lob/rn_L_RNA_scaffolder.fasta \
  -s Lobster \
  -o ./Lobster_EarlGrey \
  -t 60 \
  -i 10 \
  -f 1000 \
  -a 3 \
  -n 20 \
  -c yes \
  -m yes \
  -d yes -r arthropoda -l /worker_data/huyha/lob/lib/nr.combine.fa

# build denovo database 
pixi run BuildDatabase -name lobster_db_3 \
/worker_data/huyha/lob/assembly.5kb.fa

# Step 2: Run RepeatModeler
pixi run RepeatClassifier -consensi \
/worker_data/huyha/lob/rep/RM_310560.ThuSep251111572025/consensi.fa -stockholm \
/worker_data/huyha/lob/rep/RM_310560.ThuSep251111572025/families.stk

#Merge database from denovo and ear
cat /worker_data/huyha/lob/rep/RM_490946.ThuSep251241332025/RM_prefixed.fa \
/worker_data/huyha/lob/rep/RM_310560.ThuSep251111572025/RM_prefixed.fa \
/worker_data/huyha/lob/dnapipeTE/trim.M2_1_R1.fastq.gz/TRIN_prefixed.fa \
/worker_data/huyha/lob/dnapipeTE/trim.M2_2_R1.fastq.gz/TRIN_prefixed.fa > \
/worker_data/huyha/lob/lib
pixi run seqkit split -s 50000 -O chunks_wg  rn_L_RNA_scaffolder.fasta

# Move to working directory
cd  rp/chunks_wg || exit

# Run RepeatMasker in parallel (12 jobs × 4 threads = 48 cores total)
ls *.fasta | parallel -j 12 '
  base=$(basename {} .fasta);
  pixi run RepeatMasker \
    -pa 64 \
    -lib  lib/arthropoda_nr.combine.fa \
    -xsmall -gff -a \
    {} \
    > ${base}.log 2>&1
'
cat  rp/chunks_wg/rn_L_RNA_scaffolder.part_*.fasta.masked >  rp/Rp_element_Final/rn_L_RNA_scaffolder.masked.fasta
# Add the header to the final file
echo "repeat_family,Count,TotalLength_bp" > repeat_summary_statistics.csv
cat repeat_summary_temp.csv >> repeat_summary_statistics.csv
rm repeat_summary_temp.csv

echo "Summary file 'repeat_summary_statistics.csv' created."
cat  rp/chunks_wg/rn_L_RNA_scaffolder.part_*.fasta.out >  rp/new/rn_L_RNA_scaffolder.out
cat  rp/chunks_wg/rn_L_RNA_scaffolder.part_*.fasta.align >  rp/new/rn_L_RNA_scaffolder.align
cat  rp/chunks_wg/rn_L_RNA_scaffolder.part_*.fasta.out.gff >  rp/Rp_element_Final/rn_L_RNA_scaffolder.gff
