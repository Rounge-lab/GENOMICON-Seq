#!/bin/bash

# Check if the correct number of arguments is passed; if not, print usage
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 BAMS_DIR OUTPUT_DIR BAMS_TAGGED_DIR VCF_OUTPUT_DIR NUMBER_OF_CORES KEY_TABLE"
    exit 1
fi

# Assign arguments to variables
BAMS_DIR="$1"
OUTPUT_DIR="$2"
BAMS_TAGGED_DIR="$3"
VCF_OUTPUT_DIR="$4"
NUMBER_OF_CORES="$5"

# Define the input CSV file path
INPUT_CSV="$6"

# Ensure output directories exist
mkdir -p "${BAMS_TAGGED_DIR}" "${VCF_OUTPUT_DIR}"

# Step 1: AddOrReplaceReadGroups and index BAMs
while IFS=',' read -r sample_id original_bam_file_suffix tagged_bam_suffix sample_index copy_number mutation_rate read_number sample_number tech_rep original_sample_index normal_cancer; do
    # Remove quotes
    sample_id=$(echo $sample_id | sed 's/"//g')
    original_bam_file_suffix=$(echo $original_bam_file_suffix | sed 's/"//g')
    tagged_bam_suffix=$(echo $tagged_bam_suffix | sed 's/"//g')
    normal_cancer=$(echo $normal_cancer | sed 's/"//g')

    # Construct the full paths for the input and output BAM files
    input_bam="${BAMS_DIR}/${sample_id}${original_bam_file_suffix}"
    output_bam="${BAMS_TAGGED_DIR}/${sample_id}${tagged_bam_suffix}"

    # Add read groups
    singularity exec --bind $(pwd):/data gatk.sif gatk AddOrReplaceReadGroups \
        -I "/data/${input_bam}" \
        -O "/data/${output_bam}" \
        -RGID 1 \
        -RGLB lib1 \
        -RGPL ILLUMINA \
        -RGPU unit1 \
        -RGSM "${normal_cancer}"

    # Index the BAM file with the corrected path
    samtools index "${output_bam}"
done < <(tail -n +2 "${INPUT_CSV}") # Process the CSV file

# Step 2: Process samples by read_number and copy_number, identifying normal and cancer samples for Mutect2 analysis, with parallelization
awk -F',' '{gsub(/"/, "", $0); print $7, $5}' "${INPUT_CSV}" | sort | uniq | while read read_number copy_number; do
    normal_sample_id=$(grep ",${copy_number},.*${read_number}," "${INPUT_CSV}" | grep "normal" | cut -d',' -f1 | sed 's/"//g')
    normal_sample_suffix=$(grep ",${copy_number},.*${read_number}," "${INPUT_CSV}" | grep "normal" | cut -d',' -f3 | sed 's/"//g')
    normal_normal_cancer=$(grep ",${copy_number},.*${read_number}," "${INPUT_CSV}" | grep "normal" | cut -d',' -f11 | sed 's/"//g')

    grep ",${copy_number},.*${read_number}," "${INPUT_CSV}" | grep -v "normal" | sed 's/"//g' | while IFS=',' read -r cancer_sample_id _ cancer_sample_suffix _ _ _ _ _ _ _ cancer_normal_cancer; do
        normal_bam="${BAMS_TAGGED_DIR}/${normal_sample_id}${normal_sample_suffix}"
        cancer_bam="${BAMS_TAGGED_DIR}/${cancer_sample_id}${cancer_sample_suffix}"

        echo singularity exec --bind $(pwd):/data gatk.sif gatk Mutect2 \
            -R /data/ref_genome/chr1.fasta \
            -I "/data/${normal_bam}" \
            -I "/data/${cancer_bam}" \
            -normal "${normal_normal_cancer}" \
            --germline-resource /data/necessary_files/chr1_af_only_gnomad.hg38.vcf.gz \
            --panel-of-normals /data/necessary_files/chr1_1000g_pon.hg38.vcf.gz \
            -O "/data/${VCF_OUTPUT_DIR}/${normal_normal_cancer}_${cancer_normal_cancer}.vcf"
    done
done | parallel -j "${NUMBER_OF_CORES}"
