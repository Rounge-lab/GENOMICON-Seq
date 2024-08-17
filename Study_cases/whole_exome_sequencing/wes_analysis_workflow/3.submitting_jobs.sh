#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <path_to_config_dir> <base_slurm_script_name_without_extension> <path_to_input_data>"
    exit 1
fi

# Arguments
CONFIG_DIR="$1"            # Directory containing your configuration files
SLURM_SCRIPT_BASE="$2"     # Base name for your SLURM scripts, without the file extension
INPUT_DATA_DIR="$3"        # Path to the large input data directory

# Loop through each config file
for config_file in "$CONFIG_DIR"/config_mut_*.yml; do
    # Create a unique directory name for this job
    job_dir="$(basename "$config_file" .yml)_jobdir"
    mkdir -p "$job_dir"

    # Create symbolic links to the input data inside the job directory
    ln -s "$INPUT_DATA_DIR"/* "$job_dir"

    # Copy the SLURM script to the job directory and modify it
    temp_slurm_script="${job_dir}/${SLURM_SCRIPT_BASE}_$(basename "$config_file" .yml).slurm"
    cp "${SLURM_SCRIPT_BASE}.slurm" "$temp_slurm_script"

    # Modify the SLURM script for the new job directory
    sed -i "s|CONFIG_FILE_PLACEHOLDER|$config_file|" "$temp_slurm_script"
    sed -i "s|CURRENT_WORKING_DIR|$PWD/$job_dir|g" "$temp_slurm_script"
    sed -i "s|STORAGE_DIR_PLACEHOLDER|$PWD/$job_dir|" "$temp_slurm_script"

    # Submit the job
    sbatch "$temp_slurm_script"

    sleep 15
done

