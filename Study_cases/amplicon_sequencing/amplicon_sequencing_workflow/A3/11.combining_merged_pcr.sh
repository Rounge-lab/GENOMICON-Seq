#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_file_name>"
    exit 1
fi

input_dir=$1
output_file=$2

# Prepare the output file and write the header
echo "pcr_mutations,copy_number,pcr_mutations_table" > "${output_file}"

# Process each _merged.csv file
find "${input_dir}" -type f -name "*_merged.csv" | while read file_path; do
    table_name=$(basename "${file_path}" "_merged.csv" | sed 's/pcr_mutations_//')
    
    # Skip the header of each file and add the table name as a new column
    awk -v table_name="${table_name}" 'NR>1 {print $0","table_name}' FS=',' OFS=',' "${file_path}" >> "${output_file}"
done

echo "Merging complete. Output file: ${output_file}"

