#!/bin/bash

# Check if exactly two arguments are given
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <source_directory> <destination_directory>"
    exit 1
fi

source_dir=$1
dest_dir=$2

# Find all directories with _jobdir suffix
find "${source_dir}" -type d -name "*_jobdir" | while read jobdir; do
    # Extract jobdir name without 'config_' prefix and '_jobdir' suffix
    jobdir_name=$(basename "${jobdir}" | sed 's/config_//;s/_jobdir//')
    
    # Navigate through the required subdirectory structure
    find "${jobdir}" -type d -regex ".*/[0-9]+/Output_data_ampliseq+/sample_[0-9]+/tech_replicate_[0-9]+/PCR/PCR_[0-9]+" | while read pcr_dir; do
        # Extract sample, tech_replicate, and PCR identifiers
        sample=$(echo "${pcr_dir}" | grep -o 'sample_[0-9]\+')
        tech_replicate=$(echo "${pcr_dir}" | grep -o 'tech_replicate_[0-9]\+')
        pcr=$(basename "${pcr_dir}")
        
        # Construct the source file path
        src_file="${pcr_dir}/PCR_${pcr##*_}-polymerase_error_mutations.csv"
        
        # Construct the destination file name
        dest_file="${dest_dir}/pcr_mutations_${sample}_${jobdir_name}_${tech_replicate}_${pcr}.csv"
        
        # Copy and rename the file
        if [ -f "${src_file}" ]; then
            cp "${src_file}" "${dest_file}"
            echo "Copied to ${dest_file}"
        else
            echo "File not found: ${src_file}"
        fi
    done
done

