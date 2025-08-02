#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <master_csv_file.csv> <number_of_parts> <output_directory>"
    exit 1
fi

master_csv_file=$1
num_parts=$2
output_dir=$3


total_lines=$(wc -l < "$master_csv_file")
((lines_per_file = (total_lines + num_parts - 1) / num_parts))

awk -v lines_per_file="$lines_per_file" -v output_dir="$output_dir" '
    BEGIN { file_count = 1 }
    {
        if (NR == 1 || (NR-1) % lines_per_file == 0) {
            close(output_file);
            file_name = sprintf("%s/temp_reads_%d.csv", output_dir, file_count);
            file_count++;
            output_file = file_name;
        }
        print $0 > output_file;
    }
    END { close(output_file); }
' "$master_csv_file"

echo "Splitting complete. Files are in $output_dir"
