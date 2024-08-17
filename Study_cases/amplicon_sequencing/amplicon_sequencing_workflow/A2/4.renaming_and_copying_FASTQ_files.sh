#!/bin/bash

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <path_to_search> <genome_name> <destination_path>"
  exit 1
fi

# Starting directory is taken from the first argument
start_dir=$1
genome_name=$2
destination_path=$3  # New destination to copy the files to


mkdir -p "$destination_path"

# Start a unique identifier for appending filenames
unique_id=1

# First find all R1 files, pair them with R2, and append a unique ID
find "$start_dir" -type f -name "*_R1.fastq.gz" | while read -r R1_path; do
  R2_path=$(echo "$R1_path" | sed 's/_R1.fastq.gz/_R2.fastq.gz/')
  
  if [[ -f "$R2_path" ]]; then
    new_R1_path=$(echo "$R1_path" | sed "s/.fastq.gz/_S${unique_id}.fastq.gz/")
    new_R2_path=$(echo "$R2_path" | sed "s/.fastq.gz/_S${unique_id}.fastq.gz/")

    mv "$R1_path" "$new_R1_path"
    mv "$R2_path" "$new_R2_path"

    unique_id=$((unique_id + 1))
  fi
done

# Next, add a sample_* prefix to all files within generated_reads directories
find "$start_dir" -type f -regex ".*/sample_[0-9]+/generated_reads/.*_S[0-9]+\.fastq.gz" | while read -r file_path; do
  sample_prefix=$(echo "$file_path" | grep -o 'sample_[0-9]*')
  
  dir_path=$(dirname "$file_path")
  filename=$(basename "$file_path")
  new_file_path="${dir_path}/${sample_prefix}_${filename}"

  mv "$file_path" "$new_file_path"
done

# All the steps under change the name of the generated reads so that they will contain the info of the used paramters to make he files
# And to prepare the FASTQ file names to be compatible with the TaME-seq pipeline

# Convert all underscores to dashes in filenames, except for the one before _S* 
find "$start_dir" -type f -regex ".*/sample_[0-9]+/generated_reads/.*_S[0-9]+\.fastq.gz" | while read -r file_path; do
  dir_path=$(dirname "$file_path")
  filename=$(basename "$file_path")
  
  base_name=$(echo "$filename" | sed -r 's/(_S[0-9]+)\.fastq.gz$//')
  identifier=$(echo "$filename" | grep -o "_S[0-9]*\.fastq.gz$")
  new_base_name=$(echo "$base_name" | tr '_' '-')
  new_filename="${new_base_name}${identifier}"
  new_file_path="${dir_path}/${new_filename}"
  
  mv "$file_path" "$new_file_path"
done

# Insert genome_name into the filenames
find "$start_dir" -type f -regex ".*/sample_[0-9]+/generated_reads/.*_S[0-9]+\.fastq.gz" | while read -r file_path; do
    dir_path=$(dirname "$file_path")
    filename=$(basename "$file_path")
    
    new_filename=$(echo "$filename" | sed -r "s/(.*)(-R[12])_S([0-9]+)(\.fastq\.gz)/\1\2-${genome_name}_S\3\4/")
    
    new_file_path="${dir_path}/${new_filename}"
    
    if [ "$file_path" != "$new_file_path" ]; then
        mv "$file_path" "$new_file_path"
    fi
done

# Add -F or -R depending on PCR-1 or PCR-2
find "$start_dir" -type f -regex ".*/sample_[0-9]+/generated_reads/.*-HPV[0-9]+_S[0-9]+\.fastq.gz" | while read -r file_path; do
    dir_path=$(dirname "$file_path")
    filename=$(basename "$file_path")
    
    if [[ "$filename" =~ PCR-1 ]]; then
        new_filename=$(echo "$filename" | sed -r "s/(-HPV[0-9]+)_S([0-9]+)(\.fastq.gz)/\1-F_S\2\3/")
    elif [[ "$filename" =~ PCR-2 ]]; then
        new_filename=$(echo "$filename" | sed -r "s/(-HPV[0-9]+)_S([0-9]+)(\.fastq.gz)/\1-R_S\2\3/")
    fi

    new_file_path="${dir_path}/${new_filename}"
    
    if [ "$file_path" != "$new_file_path" ]; then
        mv "$file_path" "$new_file_path"
    fi
done

# Append _L001_R1_001 or _L001_R2_001 before .fastq.gz based on R1 or R2 in the filename
find "$start_dir" -type f -regex ".*/sample_[0-9]+/generated_reads/.*.fastq.gz" | while read -r file_path; do
    dir_path=$(dirname "$file_path")
    filename=$(basename "$file_path")

    suffix=""
    if [[ "$filename" =~ R1 ]]; then
        suffix="_L001_R1_001"
    elif [[ "$filename" =~ R2 ]]; then
        suffix="_L001_R2_001"
    fi

    new_filename="${filename%.fastq.gz}${suffix}.fastq.gz"
    new_file_path="${dir_path}/${new_filename}"

    if [ "$file_path" != "$new_file_path" ]; then
        mv "$file_path" "$new_file_path"
    fi
done

# Perform the final edit to include the modified directory portion in the filename
find "$start_dir" -type f -regex ".*/config_mut_r_.*_jobdir/.*/sample_[0-9]+/generated_reads/.*\.fastq.gz" | while read -r file_path; do
    dir_path=$(dirname "$file_path")
    filename=$(basename "$file_path")

    # Extract the part of the directory name between 'config_' and '_jobdir'
    config_part=$(echo "$dir_path" | grep -oP '(?<=config_).*(?=_jobdir)')
    # Replace underscores with dashes in the extracted part
    modified_config_part=$(echo "$config_part" | tr '_' '-')

    # Construct the new filename by inserting the modified directory part after 'sample-*'
    new_filename=$(echo "$filename" | sed -r "s/(sample-[0-9]+)-/\1-${modified_config_part}-/")
    new_file_path="${dir_path}/${new_filename}"

    if [ "$file_path" != "$new_file_path" ]; then
        mv "$file_path" "$new_file_path"
    fi
done

# Before copying files, perform modifications to remove -PCR-*-R* and change _S to -S
find "$start_dir" -type f -regex ".*/config_mut_r_.*_jobdir/.*/sample_[0-9]+/generated_reads/.*\.fastq.gz" | while read -r file_path; do
    dir_path=$(dirname "$file_path")
    filename=$(basename "$file_path")

    # Remove -PCR-*-R* pattern
    modified_filename=$(echo "$filename" | sed -r "s/-PCR-[0-9]+-R[1-2]//")

    # Change underscore before S to dash for the sequence part
    modified_filename=$(echo "$modified_filename" | sed -r "s/_(S[0-9]+_L[0-9]+_R[1-2]_001)/-\1/")

    new_file_path="${dir_path}/${modified_filename}"

    if [ "$file_path" != "$new_file_path" ]; then
        mv "$file_path" "$new_file_path"
    fi
done

# Copy all relevant files to the specified destination directory
find "$start_dir" -type f -regex ".*/config_mut_r_.*_jobdir/.*/sample_[0-9]+/generated_reads/.*\.fastq.gz" | while read -r file_path; do
    # Simply copy the file to the destination directory
    cp "$file_path" "$destination_path/"
done

echo "All files have been copied to $destination_path."
