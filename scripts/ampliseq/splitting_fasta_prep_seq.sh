#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <path_to_csv_directory> <path_to_master_fasta_file.fasta> <output_directory>"
    exit 1
fi

csv_dir=$1
master_fasta_file=$2
output_dir=$3

mkdir -p "$output_dir"

for csv_file in "$csv_dir"/temp_reads_*.csv; do
    [ -e "$csv_file" ] || continue

    awk -F ';' '{print $1}' "$csv_file" > "${csv_file%.csv}_seq_names.txt"

    file_num=$(basename "$csv_file" | cut -d '_' -f 3 | cut -d '.' -f 1)

    seqkit grep -f "${csv_file%.csv}_seq_names.txt" "$master_fasta_file" > "$output_dir/temp_$file_num.fasta"

    rm "${csv_file%.csv}_seq_names.txt"
done

rm "$master_fasta_file"

echo "FASTA files have been created in $output_dir. Master FASTA file is removed."
