#!/bin/bash

csv_list=$1
fasta_list=$2
output_fasta=$3  

while read -r file; do
    rm "$file"
done < "$csv_list"

fasta_dir=$(dirname "$(head -n 1 "$fasta_list")")
temp_fasta="$fasta_dir/temp_fasta.fasta"

while read -r file; do
    cat "$file" >> "$temp_fasta"
    rm "$file"
done < "$fasta_list"
gzip "$temp_fasta"
mv "$temp_fasta.gz" "$output_fasta"

echo "Merging and cleanup complete! Output files: $output_csv, $output_fasta"
