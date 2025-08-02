#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 file_list num_cores"
    exit 1
fi

file_list="$1"
num_cores="$2"

process_file() {
    local input_bed_gz_file="$1"
    local output_csv_file="${input_bed_gz_file%.bed.gz}_pcr.csv"

    gzip -cd "$input_bed_gz_file" | awk 'BEGIN{OFS=";"} {print $4, $2":"$3, $3-$2+1}' > "$output_csv_file"

    echo "Conversion complete. Output written to $output_csv_file"
}

export -f process_file

# Use xargs for parallel processing, specifying the number of cores
cat "$file_list" | xargs -P "$num_cores" -I {} bash -c 'process_file "{}"' 
