#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

input_file=$1
output_file=$2

awk -F ';' '{
    # Replace spaces and commas with underscores in the 5th column
    gsub(/ |,/, "_", $5);
    # Combine the 1st and 5th column with an underscore
    combined=$1"_"$5;
    # Print the combined columns and the 3rd column, separated by semicolon
    print combined " " $3
}' "$input_file" > "$output_file"

echo "Processing complete. Output saved to $output_file"
