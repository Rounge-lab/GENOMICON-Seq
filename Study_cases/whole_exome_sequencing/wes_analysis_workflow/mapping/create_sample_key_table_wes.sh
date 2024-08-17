
#!/bin/bash

# Validate input arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

input_file="$1"
tmp_output_file="tmp_$2"
output_file="$2"

# Process the input file to create a temporary output file
awk -F'\t' 'BEGIN { OFS="\t"; print "sample_id", "sample_index", "copy_number", "mutation_rate", "read_number", "sample_number", "tech_rep" }
NR > 1 {
    sample_index = substr($1, match($1, /-S[0-9]+/)+1)
    match($1, /sample-([^-]+)-mut-r-([^r]+)-r-([^-]+)-c-([^-]+)-tech-replicate-([^ -]+)-/, arr)
    print $1, sample_index, arr[4], arr[2], arr[3], arr[1], arr[5]
}' "$input_file" > "$tmp_output_file"


# Phase 1: Assign unique original_sample_index based on groups
awk -F'\t' 'BEGIN {OFS="\t"; idx=1}
{
    if(NR == 1) { print $0, "original_sample_index"; next }
    key = $3 FS $4 FS $5
    if(!(key in group)) { group[key] = idx++; }
    print $0, group[key]
}' "$tmp_output_file" > "$output_file"

# Clean up
rm "$tmp_output_file"

