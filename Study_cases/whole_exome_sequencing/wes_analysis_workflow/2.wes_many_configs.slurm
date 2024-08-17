#!/bin/bash
#SBATCH --account=ec34
#SBATCH --job-name=is_wes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=8
#SBATCH --time=3:00:00

# define the storage directory - placeholder
STORAGE_DIR="STORAGE_DIR_PLACEHOLDER"

# cleanup command to copy the entire Output_data directory from the scratch drive to the permanent storage
cleanup "cp -r $SCRATCH $STORAGE_DIR"

# Run the Singularity container with the bound necessary directories 
singularity exec \
    --bind CURRENT_WORKING_DIR:/usr/src/app/pipeline \
    --bind $SCRATCH:/usr/src/app/scratch \
    --pwd /usr/src/app/pipeline \
    genomicon-seq.sif \
    snakemake -j 1 -p -s ./Snakefile_wes --configfile CONFIG_FILE_PLACEHOLDER

