#!/bin/bash

csv_list=$1
fasta_list=$2
output_fasta=$3  

while read -r file; do
    rm "$file"
done < "$csv_list"

while read -r file; do
    cat "$file" >> temp_fasta.fasta
    rm "$file"
done < "$fasta_list"
gzip temp_fasta.fasta
mv temp_fasta.fasta.gz "$output_fasta"

echo "Cleanup complete! Removed CSV files and created output file: $output_fasta"
