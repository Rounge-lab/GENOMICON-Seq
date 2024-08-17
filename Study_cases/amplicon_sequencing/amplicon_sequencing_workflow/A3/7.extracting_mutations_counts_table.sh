#!/bin/bash

# Check for correct usage
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

input_dir="$1"
output_dir="$2"

# Ensure the output directory exists
mkdir -p "$output_dir"

# Loop through each subdirectory of the input directory
find "$input_dir" -mindepth 1 -maxdepth 1 -type d | while read sub_dir; do
    # Find the .yml file within the subdirectory
    yml_file=$(find "$sub_dir" -maxdepth 1 -type f -name "*.yml")
    if [[ ! -z "$yml_file" ]]; then
        yml_file_name=$(basename "$yml_file" .yml)
        
        # Look for numeric-named directories within the subdirectory
        find "$sub_dir" -type d -regex '.*/[0-9]+$' | while read numeric_dir; do
            # Define the path to Output_data_12 inside the numeric directory
            output_data_dir="$numeric_dir/Output_data_ampliseq"
            
            if [ -d "$output_data_dir" ]; then
                # Find all sample_* folders within Output_data_12
                find "$output_data_dir" -type d -name 'sample_*' | while read sample_dir; do
                    sample_folder_name=$(basename "$sample_dir")
                    
                    # Find the target CSV file within the sample_* folder
                    find "$sample_dir" -type f -name "*_mutation_counts.csv" | while read csv_file; do
                        csv_file_name=$(basename "$csv_file")
                        # Construct the new file name
                        new_file_name="${sample_folder_name}_${yml_file_name}_${csv_file_name}"
                        
                        # Copy and rename the file to the output directory
                        cp "$csv_file" "$output_dir/$new_file_name"
                    done
                done
            fi
        done
    fi
done
