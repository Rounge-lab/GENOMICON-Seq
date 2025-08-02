#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <csv_table> <input_fasta.gz> <output_fasta>"
    exit 1
fi

csv_table=$1
input_fasta_gz=$2
output_fasta=$3

cut -d ';' -f 1 "$csv_table" > temp_names.txt

seqkit grep --pattern-file temp_names.txt "$input_fasta_gz" > "$output_fasta"

rm "$input_fasta_gz"
rm temp_names.txt

echo "Filtering complete. Output saved in $output_fasta."
