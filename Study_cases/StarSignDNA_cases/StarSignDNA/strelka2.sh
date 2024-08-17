#!/bin/bash

# Usage: ./strelka2.sh <normal_sample.bam> <tumor_samples_dir> <reference_genome.fasta> <regions.bed> <num_cores>

NORMAL_BAM=$1              # The normal sample BAM file
TUMOR_DIR=$2               # Directory containing tumor BAM files
REFERENCE_GENOME=$3        # Reference genome file
REGIONS_FILE=$4            # BED file with restricted regions
NUM_CORES=$5               # Number of cores for parallel processing

# Iterate over each tumor BAM file in the specified directory
for TUMOR_BAM in ${TUMOR_DIR}/*.bam
do
    # Extract the base name of the file for creating unique run directories
    TUMOR_BASENAME=$(basename ${TUMOR_BAM} .bam)

    # Set up the run directory based on tumor sample name
    RUN_DIR="./strelka_workflow_${TUMOR_BASENAME}"
    mkdir -p ${RUN_DIR}

    # Configure Strelka2 somatic workflow
    configureStrelkaSomaticWorkflow.py \
        --normalBam ${NORMAL_BAM} \
        --tumorBam ${TUMOR_BAM} \
        --referenceFasta ${REFERENCE_GENOME} \
        --callRegions ${REGIONS_FILE} \
        --runDir ${RUN_DIR}

    # Run the workflow
    ${RUN_DIR}/runWorkflow.py -m local -j ${NUM_CORES}
done

