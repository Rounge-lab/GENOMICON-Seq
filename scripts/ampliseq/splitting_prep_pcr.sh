#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <file.csv> <number_of_files> <output_directory>"
    exit 1
fi

input_file="$1"
number_of_files="$2"
output_dir="$3"
base_filename=$(basename "$input_file" .csv)

header=$(head -n 1 "$input_file")

total_lines=$(($(wc -l < "$input_file") - 1))
lines_per_file=$(((total_lines + number_of_files - 1) / number_of_files))

mkdir -p "$output_dir"

current_line=2 
for ((i=1; i<=number_of_files; i++)); do
    start_line=$current_line
    end_line=$((start_line + lines_per_file - 1))

    if [ "$i" -eq "$number_of_files" ]; then
        end_line=$((total_lines + 1))
    fi

    (echo "$header"; sed -n "${start_line},${end_line}p" "$input_file") > "${output_dir}/${base_filename}_${i}.csv"

    current_line=$((end_line + 1))
done


echo "Splitting complete. Files are saved in ${output_dir}"
