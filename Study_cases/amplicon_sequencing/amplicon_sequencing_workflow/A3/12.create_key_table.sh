#!/bin/sh

# Check for the correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

# Assign arguments to variables
input_file="$1"
output_file="$2"

# Temp file for intermediate results
temp_file=$(mktemp)

# First Pass: Generate unique indexes without considering -F or -R suffix
awk 'BEGIN {
    FS=OFS="\t"
    print "sample_id", "sample_index", "mutation_rate", "copy_number", "polymerase_error_rate", "number_reads", "sample_rep", "tech_rep", "fr_suffix"
}
NR > 1 {
    # Extracting parts of the sample_id
    match($1, /sample-([^-]+)-mut-r-([e0-9.-]+)-r-([^-]+)-c-([^-]+)-poly-r-([e0-9.-]+)-tech-replicate-([^-]+)-HPV16-([FR])/, arr)
    
    # Assign extracted values
    sample_rep = arr[1]
    mutation_rate = arr[2]
    number_reads = arr[3]
    copy_number = arr[4]
    polymerase_error_rate = arr[5]
    tech_rep = arr[6]
    fr_suffix = arr[7] # Capture -F or -R suffix

    # Generate unique ID for grouping
    idWithoutHPV16 = substr($1, 1, match($1, /-HPV16/) - 1)
    if (!(idWithoutHPV16 in uniqueIds)) {
        uniqueIds[idWithoutHPV16] = length(uniqueIds) + 1
    }
    
    # Print intermediate results including fr_suffix for second pass
    print $1, uniqueIds[idWithoutHPV16], mutation_rate, copy_number, polymerase_error_rate, number_reads, sample_rep, tech_rep, fr_suffix
}' $input_file > $temp_file

# Second Pass: Append -F or -R to the unique index based on fr_suffix
awk 'BEGIN {
    FS=OFS="\t"
    print "sample_id", "sample_index", "mutation_rate", "copy_number", "polymerase_error_rate", "number_reads", "sample_rep", "tech_rep"
}
NR > 1 {
    # Reconstruct the sample_index with -F or -R suffix
    $2 = "S_" $2 "-" $9
    
    # Print the final line without the fr_suffix column
    print $1, $2, $3, $4, $5, $6, $7, $8
}' $temp_file > $output_file

# Clean up
rm -f $temp_file

echo "Processing complete. Output saved to $output_file"

