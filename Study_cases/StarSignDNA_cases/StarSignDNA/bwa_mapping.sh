#!/bin/bash

# Check if sufficient arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <reference path> <fastq dir> <output dir>"
    exit 1
fi

# Assign command-line arguments to variables
REF="$1"
DIR="$2"
OUTDIR="$3"


# Loop through all R1 files
for R1 in ${DIR}*_R1.fastq.gz; do
    # Corresponding R2 file
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
    
    # Extract the sample base name by removing the directory prefix and suffix
    SAMPLE=$(basename ${R1} _R1.fastq.gz)

    # Run BWA MEM
    bwa mem -t 4 $REF $R1 $R2 | samtools view -bS - > ${OUTDIR}${SAMPLE}.bam

    # Optionally, you can sort the BAM files immediately
    samtools sort -o ${OUTDIR}${SAMPLE}_sorted.bam ${OUTDIR}${SAMPLE}.bam

    # And index the sorted BAM
    samtools index ${OUTDIR}${SAMPLE}_sorted.bam
done

echo "BWA processing complete."

