#!/bin/bash
source config.sh

echo "Starting 7. Finding candidate sex-specific and positive control marker"
# Exit immediately if any command fails
set -e 

# --- Configuration ---
BASE_DIR="/mnt/10T/lobster_project"
KMC_DIR="$BASE_DIR/kmer/kmer_c"
TMP_DIR="$KMC_DIR/tmp"
SCRIPT_DIR="$BASE_DIR/script_project"
FINAL_OUT_DIR="$BASE_DIR/kmer"

# Create directories
mkdir -p "$TMP_DIR"
mkdir -p "$FINAL_OUT_DIR"

# Define k-mer sizes
kmer_sizes=(21 33 55)

# --- 1. K-mer Counting ---
echo "Starting KMC counting..."

run_kmc() {
    input_dir=$1
    prefix=$2
    
    # Check if directory exists
    if [ ! -d "$input_dir" ]; then
        echo "Error: Directory $input_dir does not exist."
        return
    fi

    # Filter only FASTQ files (gz or uncompressed) to avoid reading junk files
    find "$input_dir" -maxdepth 1 -name "*.fq.gz" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fastq" | while read i; do
        o=$(basename "$i")
        # Cleanup filename
        o=${o/merge_/}
        o=${o/trim./}
        o=${o/.fq.gz/}
        o=${o/.fastq.gz/}

        for k in "${kmer_sizes[@]}"; do
            output_db="$KMC_DIR/kmc.$prefix$o.k$k"
            
            # Skip if already exists to save time (Optional)
            if [ ! -d "$output_db" ]; then
                echo "Processing $o with k=$k..."
                kmc -k"$k" -t32 -m80 -ci2 "$i" "$output_db" "$TMP_DIR"
            else
                echo "Skipping $o (k=$k), output already exists."
            fi
        done
    done
}

# Run KMC
run_kmc "$BASE_DIR/data/merge" ""
run_kmc "$BASE_DIR/DarT-seq_data/trim" "dart."


# --- 2. Merging K-mers ---
# STOP: Ensure your .operate files inside $SCRIPT_DIR point to the correct files generated above!
echo "Merging k-mers..."

# You might need to generate these operation files dynamically if filenames vary.
# Assuming your manually created files are correct:
kmc_tools complex "$SCRIPT_DIR/33.kmer_operate.F"
kmc_tools complex "$SCRIPT_DIR/33.kmer_operate.M"
kmc_tools complex "$SCRIPT_DIR/55.kmer_operate.F"
kmc_tools complex "$SCRIPT_DIR/55.kmer_operate.M"


# --- 3. Extract Unique K-mers (Subtraction) ---
echo "Extracting unique k-mers..."

# Function to subtract, dump, and convert to FASTA
process_unique() {
    local k=$1
    local set1=$2 # Target (e.g., F)
    local set2=$3 # Subtract (e.g., M)
    
    local in_db1="$KMC_DIR/${set1}.merged_kmc.k$k"
    local in_db2="$KMC_DIR/${set2}.merged_kmc.k$k"
    local uniq_db="$KMC_DIR/uniq_${set1}_kmc.k$k"
    local dump_txt="$KMC_DIR/3.uniq_${set1}_kmc.k$k.txt"
    local out_fasta="$FINAL_OUT_DIR/k$k.${set1}.merged_kmers.fasta"

    # Subtract
    kmc_tools simple "$in_db1" "$in_db2" kmers_subtract "$uniq_db" -ci3
    
    # Transform to text
    kmc_tools transform "$uniq_db" dump "$dump_txt" -ci3
#WGS
mkdir -p /mnt/10T/lobster_project/kmer/clean_read
mkdir -p /mnt/10T/lobster_project/kmer/clean_read/QC

for i in `ls /mnt/10T/lobster_project/kmer/clean_read/k33.clean_clean_merge_F*_R1.fastq.gz`;do
a=`echo $i | sed -e  's/R1/R2/'`
b=`basename ${a}`
c=`basename ${i}`
bbduk.sh in1=$i in2=$a ref=/mnt/10T/lobster_project/kmer/k55.M.merged_kmers.fasta -Xmx90g  usejni=t \
         outm1=/mnt/10T/lobster_project/kmer/clean_read/k55_$c outm2=/mnt/10T/lobster_project/kmer/clean_read/k55_$b k=21 hdist=0 \
         stats=$b.stats.txt minlength=50 &>  k55.M.$c.bbduck.log
done

for i in `ls /mnt/10T/lobster_project/kmer/clean_read/k33.clean_clean_merge_M*_R1.fastq.gz`;do
a=`echo $i | sed -e  's/R1/R2/'`
b=`basename ${a}`
c=`basename ${i}`
bbduk.sh in1=$i in2=$a ref=/mnt/10T/lobster_project/kmer/k55.F.merged_kmers.fasta -Xmx90g  usejni=t\
         outm1=/mnt/10T/lobster_project/kmer/clean_read/k55_$c outm2=/mnt/10T/lobster_project/kmer/clean_read/k55_$b k=21 hdist=0 \
         stats=$b.stats.txt minlength=50 &>  k55.F.$c.bbduck.log
done

# DarT-seq
mkdir /mnt/10T/lobster_project/kmer/DarT_seq_clean_read
rm -f bbduk.F.txt 
for i in /mnt/10T/lobster_project/kmer/DarT_seq_clean_read/k33.clean_trim.Female*.fq.gz; do
    c=$(basename "$i")  
    echo $i
    [[ -s "/mnt/10T/lobster_project/kmer/DarT_seq_clean_read/k55_$c" ]] || \
    bbduk.sh in="$i" ref=/mnt/10T/lobster_project/kmer/k55.F.merged_kmers.fasta -Xmx75g usejni=f \
             outm=/mnt/10T/lobster_project/kmer/DarT_seq_clean_read/k55_"$c" \
             k=21 hdist=0 stats="$c.k55.stats.txt" minlength=50 threads=8 &> k55.M."$c".bbduk.log >> bbduk.F.txt
done


rm -f bbduk.M.txt 
for i in  /mnt/10T/lobster_project/kmer/DarT_seq_clean_read/k33.clean_trim.Male*.fq.gz; do
    c=`basename ${i}`   
    echo $i
    [[ -s "/mnt/10T/lobster_project/kmer/DarT_seq_clean_read/k55_$c" ]] || \
    bbduk.sh in=$i ref=/mnt/10T/lobster_project/kmer/k55.M.merged_kmers.fasta -Xmx75g usejni=f \
             outm=/mnt/10T/lobster_project/kmer/DarT_seq_clean_read/k55_"$c" \
             k=21 hdist=0 stats=$c.k55.stats.txt minlength=50 threads=8  &> k55.M.$c.bbduck.log >> bbduk.M.txt
done


parallel -j 2 "pigz -dc -p 16 /mnt/10T/lobster_project/kmer/clean_read/k55_k33.clean_clean_merge_F*_{}.fastq.gz| pigz -p 8 > /mnt/10T/lobster_project/kmer/merged_read_filtered/Merged_F.{}.fastq.gz
" ::: R1 R2 
parallel -j 2 "pigz -dc -p 16 /mnt/10T/lobster_project/kmer/clean_read/k55_k33.clean_clean_merge_M*_{}.fastq.gz| pigz -p 8 > /mnt/10T/lobster_project/kmer/merged_read_filtered/Merged_M.{}.fastq.gz
" ::: R1 R2 

mkdir -p  kmer/megahit_output/

# assembly FeMale
megahit \
    -1  kmer/merged_read_filtered/Merged_F.R1.fastq.gz \
    -2  kmer/merged_read_filtered/Merged_F.R2.fastq.gz \
    --k-min 27 --k-max 141 --k-step 20 \
    --min-contig-len 300 \
    --num-cpu-threads 40 \
    --memory 0.95 --mem-flag 0 \
    --continue \
    -o  kmer/megahit_output/F 2>  kmer/megahit_output/F.error 

exit

# assembly Male
megahit \
    -1  kmer/merged_read_filtered/Merged_M.R1.fastq.gz \
    -2  kmer/merged_read_filtered/Merged_M.R2.fastq.gz \
    --k-min 27 --k-max 141 --k-step 20 \
    --min-contig-len 300 \
    --num-cpu-threads 40\ \
    --memory 0.95 \
    --continue \
    -o  kmer/megahit_output/M 2>  kmer/megahit_output/M.error

# Rename and index reference genome
mv /mnt/10T/lobster_project/kmer/megahit_output/F/final.contigs.fa /mnt/10T/lobster_project/kmer/megahit_output/F/FEMALE_ref.fasta
bwa index /mnt/10T/lobster_project/kmer/megahit_output/F/FEMALE_ref.fasta

echo Step1
maping
parallel -j 3 '
    i={}
    a=$(echo "$i" | sed "s/_R1.fastq.gz/_R2.fastq.gz/")
    b=$(basename "$i" | sed "s/k55_k33.clean_clean_merge_//" | sed "s/_R1.fastq.gz//")
    bwa mem -t 40 -M -T 50 -B 5 -O 10 -E 3 -Y "/mnt/10T/lobster_project/kmer/megahit_output/F/FEMALE_ref.fasta" "$i" "$a" | \
    samtools view -bS ->  "/mnt/10T/lobster_project/kmer/bwa_result/WGS/$b.F.bam"
' :::  /mnt/10T/lobster_project/kmer/clean_read/k55_k33.clean_clean_merge_F3_R1.fastq.gz /mnt/10T/lobster_project/kmer/clean_read/k55_k33.clean_clean_merge_M1_R1.fastq.gz /mnt/10T/lobster_project/kmer/clean_read/k55_k33.clean_clean_merge_M2_R1.fastq.gz




echo Step2
parallel -j 3 '
    i={}
    b=$(basename "$i" | sed "s/k55_k33.clean_clean_merge_//" | sed "s/_R1.fastq.gz//")
    samtools sort -@ 27 "/mnt/10T/lobster_project/kmer/bwa_result/WGS/$b.F.bam" -o "/mnt/10T/lobster_project/kmer/bwa_result/WGS/sort.${b}.F.bam"
' ::: /mnt/10T/lobster_project/kmer/clean_read/k55_k33.clean_clean_merge_F3_R1.fastq.gz /mnt/10T/lobster_project/kmer/clean_read/k55_k33.clean_clean_merge_M1_R1.fastq.gz /mnt/10T/lobster_project/kmer/clean_read/k55_k33.clean_clean_merge_M2_R1.fastq.gz


echo Step3
#index.sort
parallel -j 6 '
    i={}
    a=$(echo "$i" | sed "s/_R1.fastq.gz/_R2.fastq.gz/")
    b=$(basename "$i" | sed "s/k55_k33.clean_clean_merge_//" | sed "s/_R1.fastq.gz//")

    samtools index "/mnt/10T/lobster_project/kmer/bwa_result/WGS/sort.${b}.F.bam"
' ::: /mnt/10T/lobster_project/kmer/clean_read/k55_k33.clean_clean_merge_F3_R1.fastq.gz /mnt/10T/lobster_project/kmer/clean_read/k55_k33.clean_clean_merge_M1_R1.fastq.gz /mnt/10T/lobster_project/kmer/clean_read/k55_k33.clean_clean_merge_M2_R1.fastq.gz



echo Step4
#maping
parallel -j 4 '
    b=$(basename {} | sed "s/k33.clean_trim.//" | sed "s/.fq.gz//")
    bwa mem -t 20 -M -T 50 -B 5 -O 10 -E 3 -Y "/mnt/10T/lobster_project/kmer/megahit_output/F/FEMALE_ref.fasta" {} | \
    samtools view -bS -> "/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/${b}.F.bam"
    
' ::: /mnt/10T/lobster_project/kmer/DarT_seq_clean_read/k33.*.fq.gz
echo Step5
#sort
parallel -j 4 '
    b=$(basename {} | sed "s/k33.clean_trim.//" | sed "s/.fq.gz//")
    samtools sort -@ 20 "/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/${b}.F.bam" -o "/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.${b}.F.bam"
    
' ::: /mnt/10T/lobster_project/kmer/DarT_seq_clean_read/k33.*.fq.gz
echo Step6
#index.sort
parallel -j 8 '
    b=$(basename {} | sed "s/k33.clean_trim.//" | sed "s/.fq.gz//")
    samtools index "/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.${b}.F.bam"
' ::: /mnt/10T/lobster_project/kmer/DarT_seq_clean_read/k33.*.fq.gz

echo Step7
samtools depth -m 100000 -aa \
/mnt/10T/lobster_project/kmer/bwa_result/WGS/sort.F1.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/WGS/sort.F2.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/WGS/sort.F3.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_01_1.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_01_2.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_02_1.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_02_2.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_03_1.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_03_2.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_04_1.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_04_2.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_05_1.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_05_2.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_06_1.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_06_2.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_07_1.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_07_2.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_08_1.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_08_2.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_09_1.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_09_2.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_10.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_11.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_12.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Female_15.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/WGS/sort.M1.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/WGS/sort.M2.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/WGS/sort.M3.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Male_01.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Male_02.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Male_03.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Male_04.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Male_05.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Male_06.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Male_07.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Male_08.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Male_09.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Male_10.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Male_11.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Male_12.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Male_14.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Male_15_1.F.bam \
/mnt/10T/lobster_project/kmer/bwa_result/dart-seq/sort.Male_15_2.F.bam \
> /mnt/10T/lobster_project/kmer/bwa_result/FEMALE.depth


echo Step8
# Note: ssp2/step2.pl is cloned from the ssp2 repository by fengtong-bio.
# Please cite: https://github.com/fengtong-bio/ssp2
perl /mnt/10T/lobster_project/canu_kmer/bwa/ssp2/step2.pl /mnt/10T/lobster_project/kmer/megahit_output/F/FEMALE_ref.fasta /mnt/10T/lobster_project/kmer/megahit_output/F/MALE_ref.fasta 25 18 /mnt/10T/lobster_project/kmer/bwa_result/FEMALE.depth /mnt/10T/lobster_project/kmer/bwa_result/MALE.depth 0.9 20 100


#Finding positive control 
awk \
-v G1_NUM=25 \
-v G2_NUM=18 \
-v GAP_LENGTH=10 \
-v MIN_LENGTH=10 \
-f  PS_region/script.awk \
 kmer/bwa_result/FEMALE.depth 

