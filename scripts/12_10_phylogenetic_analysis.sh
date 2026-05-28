#!/bin/bash
source config.sh

echo "Starting 10. Phylogenetic analysis"
grep -v "#" BUSCO/run_arthropoda_odb12/full_table.tsv | awk '$2=="Complete" {print $1}' > treefam_out/buscoID.list

treefam_scan=/path/to/treefam_tools/v1/treefam_scan/treefam_scan.pl
hmmlib=/path/to/treefam/hmm_lib

for id in `cat treefam_out/buscoID.list`; do
  fa=BUSCO/run_arthropoda_odb12/busco_sequences/single_copy_busco_sequences/${id}.faa
  perl $treefam_scan -fasta $fa -dir $hmmlib -hmm_file treefam9.hmm3 -cpu 24 -outfile treefam_out/${id}.treefamscan.tbl
done
mkdir -p phylogenicTree/fasta phylogenicTree/tbl phylogenicTree/4phylo

# Extract T.aus ortholog proteins
for i in `xlsx2csv treefam_out/Treefam_Ortho_Genes.xlsx | awk -F"," 'NR>1 {print $6}'`; do
	fa=`find BUSCO/run_arthropoda_odb12/busco_sequences/single_copy_busco_sequences -name "${i}.faa"`
	cat $fa >> phylogenicTree/fasta/T.aus_ortho.faa
done

# Extract target orthologous sequences across reference crustacean species
for i in `awk 'NR>1 {print $2}' crustacea_ref/Refseq/Crustacea.tsv | sort | uniq`; do
    blastp -query phylogenicTree/fasta/T.aus_ortho.faa -subject crustacea_ref/Refseq/ncbi_dataset/data/$i/protein.faa -outfmt 6 -max_target_seqs 1 -out phylogenicTree/tbl/$i.tsv 
    cut -f2 phylogenicTree/tbl/$i.tsv | sort | uniq > phylogenicTree/tbl/$i.list
    samtools faidx crustacea_ref/Refseq/ncbi_dataset/data/$i/protein.faa -r phylogenicTree/tbl/$i.list > phylogenicTree/fasta/$i.faa
done

# Concatenate proteins and merge
for i in `awk 'NR>1 {print $2}' crustacea_ref/Refseq/Crustacea.tsv | sort | uniq`; do
	cat phylogenicTree/fasta/$i.faa | sed -e '1!{/^>.*/d;}' | sed ':a;N;$!ba;s/\n//2g' | sed '1!s/.\{80\}/&\n/g' > phylogenicTree/fasta/$i.concanated.faa
	b=`grep ">" phylogenicTree/fasta/$i.concanated.faa`
	sed -i -e "s/$b/>$i/g" phylogenicTree/fasta/$i.concanated.faa
done

cat phylogenicTree/fasta/T.aus_ortho.faa | sed -e '1!{/^>.*/d;}' | sed ':a;N;$!ba;s/\n//2g' | sed '1!s/.\{80\}/&\n/g' > phylogenicTree/fasta/T.aus_ortho.concanated.faa
b=`grep ">" phylogenicTree/fasta/T.aus_ortho.concanated.faa`
sed -i -e "s/$b/>T.aus/g" phylogenicTree/fasta/T.aus_ortho.concanated.faa
sed -i -e "s/*//g" phylogenicTree/fasta/T.aus_ortho.concanated.faa

cat phylogenicTree/fasta/*.concanated.faa > phylogenicTree/4phylo/merged.faa
# Alignment with MAFFT
mafft --thread 48 --auto --clustalout --reorder phylogenicTree/4phylo/merged.faa > phylogenicTree/4phylo/merged.aln

# Construct phylogenetic tree using IQ-TREE
iqtree -s phylogenicTree/4phylo/merged.aln -nt 48 -m Blosum62
iqtree -s phylogenicTree/4phylo/merged.aln -b 100 -nt 48 -m Blosum62
