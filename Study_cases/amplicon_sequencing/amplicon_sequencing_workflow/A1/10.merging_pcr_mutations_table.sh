#!/bin/bash

# Check if the directory is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <directory_with_csv_files>"
    exit 1
fi

directory=$1

# Process each file
find "${directory}" -type f -name "pcr_mutations*.csv" | sed 's/_PCR_[12].csv$//' | sort | uniq | while read base_name; do
    file1="${base_name}_PCR_1.csv"
    file2="${base_name}_PCR_2.csv"
    merged_file="${base_name}_merged.csv"

    if [[ -f "$file1" && -f "$file2" ]]; then
        echo "Merging ${file1} and ${file2} into ${merged_file}"

        # Merge files and sum copy_number for identical pcr_mutations
        awk -F, 'NR==1{if(!seen){print; seen=1}} NR>1{data[$1]+=$2} END{for (mutation in data) print mutation "," data[mutation]}' "${file1}" "${file2}" > "${merged_file}"
    fi
done

