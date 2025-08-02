#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 file_list exome_bed ncores matching_length"
    exit 1
fi

file_list="$1"
exome_bed="$2"
ncores="$3"
matching_length="$4"

process_file() {
    local bed_file="$1"
    local exome_bed="$2"
    local matching_length="$3"
    local output_filtered_bed="${bed_file%.bed.gz}_filtered.bed.gz"

    bedtools intersect -a "$bed_file" -b "$exome_bed" -wo | \
        awk -v threshold="$matching_length" '$NF >= threshold {print $1, $2, $3, $4}' | \
        sort | uniq | gzip > "$output_filtered_bed"

    echo "Filtering complete for $output_filtered_bed"

    rm "$bed_file"
    echo "Removed $bed_file"
}

export -f process_file
export exome_bed
export matching_length

cat "$file_list" | xargs -P "$ncores" -I {} bash -c 'process_file "{}" "$exome_bed" "$matching_length"'
