#!/bin/bash

# Check if an input directory is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_input_folder>"
    exit 1
fi

input_dir=$(realpath "$1")
# Assuming the ref_genome directory is at the same level as the script
ref_genome_dir=$(realpath "$(dirname "$0")/ref_genome")

# Loop over each VCF file in the input directory
for vcf_file in "$input_dir"/*.vcf; do
    # Check if file exists since glob might not match anything
    if [ ! -f "$vcf_file" ]; then
        echo "No VCF files found in $input_dir."
        exit 1
    fi

    # Extract the base name of the file without the extension
    base_name=$(basename "$vcf_file" .vcf)

    # Filter the VCF file
    singularity exec --bind "$input_dir":/data --bind "$ref_genome_dir":/ref_genome gatk.sif gatk FilterMutectCalls \
        -V /data/"$base_name.vcf" \
        -O /data/"${base_name}_filtered.vcf" \
        -R /ref_genome/chr1.fasta

    # Convert the filtered VCF file to a CSV file
    singularity exec --bind "$input_dir":/data --bind "$ref_genome_dir":/ref_genome gatk.sif gatk VariantsToTable \
        -V /data/"${base_name}_filtered.vcf" \
        -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER \
        -GF GT -GF AD -GF DP -GF AF \
        --show-filtered true \
        -O /data/"${base_name}".csv
done

echo "Processing completed."

