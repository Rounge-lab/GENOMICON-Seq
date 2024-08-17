#!/bin/bash
#SBATCH --account=ec34
#SBATCH --job-name=is_ampliseq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=116
#SBATCH --mem-per-cpu=2G
#SBATCH --time=2:00:00


# placeholder the storage directory
STORAGE_DIR="STORAGE_DIR_PLACEHOLDER"

# cleanup command to copy the entire Output_data directory
# from the scratch drive to permanent storage
cleanup "cp -r $SCRATCH $STORAGE_DIR"

# Run the Singularity container with bound folders
singularity exec \
    --bind CURRENT_WORKING_DIR:/usr/src/app/pipeline \
    --bind $SCRATCH:/usr/src/app/scratch \
    --pwd /usr/src/app/pipeline \
    genomicon-seq.sif \
    snakemake -j 1 -p -s ./Snakefile_isampliseq --configfile CONFIG_FILE_PLACEHOLDER

