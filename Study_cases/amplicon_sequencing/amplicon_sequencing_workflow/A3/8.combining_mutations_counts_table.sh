#!/bin/bash

# Check if at least one argument was provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <directory_path> [output_file_name]"
    exit 1
fi

directory=$1
output_file="${2:-combined_csv.csv}" # Use the second argument as output file name or default to combined_csv.csv
data_file="data_temp.csv"
header_file="header_temp.csv"
> "$data_file" # Create or clear the temporary data file

# Extract and prepare the header from the first CSV file
head -1 "$(find "$directory" -type f -name '*.csv' | head -1)" | awk 'BEGIN{FS=OFS=","}{print "csv_table," $0}' > "$header_file"

# Process each CSV file in the directory
for csv_file in "$directory"/*.csv; do
    # Extract filename without path and extension for the new column
    csv_table_name=$(basename "$csv_file" .csv)

    # Skip the header and add the "csv_table" column
    awk -v name="$csv_table_name" 'BEGIN{FS=OFS=","} NR>1{print name, $0}' "$csv_file" >> "$data_file"
done

# Combine header and data into the final output file
cat "$header_file" "$data_file" > "$output_file"

# Clean up temporary files
rm -f "$data_file" "$header_file"

echo "Combined CSV created as $output_file"
