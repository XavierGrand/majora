#!/bin/bash

# Check that all arguments all correctly passed to the command
if [ "$#" -ne 5 ]; then
    echo "Five arguments are needed, you gave $# arguments" 
    echo "Arguments should be in the following order: reference.fasta reference_id.txt calling.vcf PRIMER_FWD PRIMER_REV"
    exit 1
else
    echo "All arguments passed correctly"
fi


# === PARAMETERS ===
FASTA="$1"  # barcodetoMap.fasta 
ID_FILE="$2" # barcode_best_ref.txt      
VCF_FILE="$3" # barcode_filtered.vcf
OUT="alignment_positions.txt"         

PRIMER_FWD="$4"
PRIMER_REV="$5"


# === Blast primers on reference sequence to get alignment positions ===

# Create database
makeblastdb -in "$FASTA" -dbtype nucl -out "$FASTA".db

# Blast forward primer
echo ">query_sequence" > query_FWD_primer.fasta
echo "$PRIMER_FWD" >> query_FWD_primer.fasta
FWD=$(blastn -query query_FWD_primer.fasta -db "$FASTA".db -evalue 1e-3 -word_size 4 -outfmt "6 sstart send" -max_target_seqs 1)
echo "forward $FWD" > alignment_res_forward.txt
# Only keep start position of the first alignment location 
cat alignment_res_forward.txt | paste -s -d '\t' | awk '{print $2 -1}' > positions_for_fasta.txt

# Blast reverse primer
echo ">query_sequence" > query_REV_primer.fasta
echo "$PRIMER_REV" >> query_REV_primer.fasta
REV=$(blastn -query query_REV_primer.fasta -db "$FASTA".db -evalue 1e-3 -word_size 4 -outfmt "6 sstart send" -max_target_seqs 1)
echo "reverse $REV" > alignment_res_reverse.txt
# Only keep start position of the second alignment location 
cat alignment_res_reverse.txt | paste -s -d '\t' | awk '{print $4}' >> positions_for_fasta.txt
# Switch file disposition from 3 rows to 1 unique row with 3 tab separated columns
paste -s -d '\t' positions_for_fasta.txt > primer_positions.txt


# === Store positions as variables ===

# Reference sequence 's ID'
ID_REF=$(awk 'NR==1 {print "gnl|hbvnuc|"$1}' $ID_FILE)

# Count nucletotides containted in concatenated reference sequence
REF_SEQ_LEN=$(cat $FASTA | sed '1d' | paste -s -d '' | wc -m)

# Calculate actual reference sequence length as: concatenated_ref_seq_length/2
REF_SEQ_LEN=$(((REF_SEQ_LEN - 1) / 2))
echo "The reference sequence $ID_REF length is: $REF_SEQ_LEN bp." >> summary.txt

START_FWD_PRIMER=$(awk '{print $1}' primer_positions.txt)
START_REV_PRIMER=$(awk '{print $2}' primer_positions.txt)
echo "The starting position of the forward primer is $START_FWD_PRIMER and of reverse primer is: $START_REV_PRIMER." >> summary.txt
# Calculate gap length resulting from PCR as: reference_sequence_length - pcr_product_length
GAP=$((REF_SEQ_LEN - (START_REV_PRIMER - START_FWD_PRIMER)))
echo "The gap resulting from PCR design as a length of: $GAP bp." >> summary.txt


# === Re-indexing vcf file according to primer alignment === 

# Considering first A of the EcoRI cutting site of HBV genome as the first base of the reference sequence
    # 1. HBV_db sequences are stored with first T or C, depending on EcoRI sequence, considered as first base. Then in concatenated reference sequence, the legnth of the actual reference sequence gives the position of the first T(or C). Thus, position of the first A is given by the following: REFERENCE_SEQUENCE_LENGTH - 2
FST_NT=$(($REF_SEQ_LEN -2))
    # 2. To set the new location of DNA
AFTER_GAP=$((START_REV_PRIMER - $FST_NT + $GAP))

echo "FST_NT: $FST_NT"
echo "START_REV_PRIMER: $START_REV_PRIMER"
echo "GAP: $GAP"

    # 3. To reposition second part of the sequence
awk -v fst_nt="$FST_NT" -v start_rev="$START_REV_PRIMER" -v gap="$GAP" -v start_fwd="$START_FWD_PRIMER" '!/^#/ {if ($2 < fst_nt) {$2 = $2} else {$2 = $2 - fst_nt}}  {print} ' $VCF_FILE > new_vcf.vcf

grep "^#" new_vcf.vcf > header.vcf
grep -v "^#" new_vcf.vcf | sort -k2,2n > new_sorted.vcf
sed -i 's/ /\t/g' new_sorted.vcf
cat header.vcf new_sorted.vcf > sorted_reindexed.vcf
rm header.vcf new_sorted.vcf


# == Extract SNP/INDEL positions, reference allele and alternative allele for visualization
awk 'BEGIN {OFS="\t"} !/^#/ {print $2, $4, $5}' sorted_reindexed.vcf > vcf.tsv


# == Add variation type to the tsv ==
awk 'BEGIN  {OFS="\t"} {A=(length($2) - length($3))} {if (A>0) {print $1, $2, $3, "deletion"} else if (A==0) {print $1, $2, $3, "SNP"} else if (A<0) {print $1, $2, $3, "insertion"}}' vcf.tsv > 02_vcf.tsv

sed -i '1iPosition\tRef\tAlt\tVariation' 02_vcf.tsv

# == Extract depth values for each variant

    # Starting with extraction of the number of ambiguous reads
awk '!/^#./ {print $8}' $VCF_FILE | awk -F";" '{print $1, $6}' | awk -F"=" '{print $2, $3}' | awk -F" " '{print $1, $3}' | awk -F" " '{print $1}' | awk -F"," '{amb= $1+$2} {print amb}' >> amb_values.txt

    # now extraction of the number of reads matching with the reference allele
awk '!/^#./ {print $8}' $VCF_FILE | awk -F";" '{print $1, $6}' | awk -F"=" '{print $2, $3}' | awk -F" " '{print $1, $3}' | awk -F" " '{print $2}' | awk -F"," '{ref_reads=$1+$2} {print ref_reads}' >> ref_reads.txt

    # now extraction of the number of reads matching with the alternative allele
awk '!/^#./ {print $8}' $VCF_FILE | awk -F";" '{print $1, $6}' | awk -F"=" '{print $2, $3}' | awk -F" " '{print $1, $3}' | awk -F" " '{print $2}' | awk -F"," '{alt_reads=$3+$4} {print alt_reads}' >> alt_reads.txt

paste amb_values.txt ref_reads.txt alt_reads.txt > oskour.csv

awk '{sum=$1+$2+$3} {print ($1/sum)*100, ($2/sum)*100, ($3/sum)*100}' oskour.csv > aled.csv

sed -i 's/-nan/none/g' aled.csv
sed -i 's/ /\t/g' aled.csv
sed -i '1iAmbiguous_reads\tRef_reads\tAlt_reads' aled.csv

paste 02_vcf.tsv aled.csv > final_vcf.tsv

# == Extract gap positions for visualization ==
START_NEW_GAP=$(($START_REV_PRIMER - $FST_NT))
END_NEW_GAP=$(($START_NEW_GAP + $GAP))
echo -e "Start_gap\t$START_NEW_GAP" > gap_positions.tsv
echo -e "End_gap\t$END_NEW_GAP" >> gap_positions.tsv
echo -e "The gap is located between $START_NEW_GAP bp and $END_NEW_GAP bp." >> summary.txt

# == Creating corrected fasta file from the concatenated reference
printf "$ID_REF" | paste -d'\t' - primer_positions.txt > ref_seq_primer2primer.bed
bedtools getfasta -fi "$FASTA" -bed ref_seq_primer2primer.bed -fo ref_seq_primer2primer.fasta

rm -f query_FWD_primer.fasta alignment_res_forward.txt query_REV_primer.fasta alignment_res_reverse.txt positions_for_fasta.txt primer_positions.txt new_vcf.vcf sorted_reindexed.vcf vcf.tsv amb_values.txt ref_reads.txt alt_reads.txt oskour.csv aled.csv *.fasta.db* *.fasta.fai 02_vcf.tsv