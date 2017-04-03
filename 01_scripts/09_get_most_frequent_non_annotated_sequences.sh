#!/bin/bash
# For each sample, eliminate the sequences with annotation found with vsearch,
# find the 100 most frequent sequences and put them in a fasta file for further
# blasts on NCBI nr/nt.

# Parameters
NUM_NON_ANNOTATED_SEQ=$1
CHIMERAFOLDER="08_chimeras"

# Get names of unwanted sequences from 09_vsearch
ls -1 09_vsearch/ | cut -d "_" -f 1 | sort -u | while read i
do
    cat 09_vsearch/"$i"_* |
        awk '{print $1}' |
        cut -d ";" -f 1 > 11_non_annotated/"$i"_with_result.ids

    # Sort them by decreasing order of count (most frequent sequences first)
    cp "$CHIMERAFOLDER"/"$i"_*_unique.fasta 11_non_annotated/"$i"_temp.fasta
    ./01_scripts/util/fasta_sort_by_count.py \
        11_non_annotated/"$i"_temp.fasta \
        11_non_annotated/"$i"_unique.fasta
    rm 11_non_annotated/"$i"_temp.fasta
done

# Remove these from "$CHIMERAFOLDER"/*.fasta and output a fasta file per sample
for i in 11_non_annotated/*.ids
do
    idfile=$(basename "$i")
    input="${idfile%_with_result.ids}"
    output="${idfile%_with_result.ids}"_without_result.fasta.gz
    ./01_scripts/util/fasta_remove.py 11_non_annotated/"$input"_unique.fasta \
        11_non_annotated/"$input"_with_result.ids \
        11_non_annotated/"$input"_without_result.fasta
done

# Get the 100 most represented unique sequences per sample
head -n $[ $NUM_NON_ANNOTATED_SEQ * 2 ] 11_non_annotated/*_without_result.fasta |
    grep -v "^==" |
    grep -vE "^$" > 12_results/most_frequent_non_annotated_sequences.fasta
