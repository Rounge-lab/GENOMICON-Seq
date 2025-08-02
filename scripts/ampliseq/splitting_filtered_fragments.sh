#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file.csv> <number_of_parts>"
    exit 1
fi

input_file="$1"
num_parts="$2"

total_lines=$(wc -l < "$input_file")

((lines_per_file = (total_lines + num_parts - 1) / num_parts))

dir=$(dirname "$input_file")

awk -v lines_per_file="$lines_per_file" -v dir="$dir" -v num_parts="$num_parts" '
    BEGIN { file_count = 1 }
    (NR-1) % lines_per_file == 0 {
        if (file_count <= num_parts) {
            output_file = dir"/worker_"file_count"_temp.csv";
            file_count++;
        }
    }
    { print $0 > output_file }
' "$input_file"

rm "$input_file"

echo "Splitting complete. Input file removed. Output files are in $dir"
