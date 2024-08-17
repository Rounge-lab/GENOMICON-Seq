#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <index_table> <table_to_process>"
    exit 1
fi

index_table=$1
table_to_process=$2

# Check if the table to be processed has a header with 'sample_id'
header=$(head -n 1 $table_to_process | grep -o 'sample_id')

if [ -n "$header" ]; then
    # Scenario 1: The table has a header, including 'sample_id'
    
    # Generate an AWK script to perform the replacement
    awk_script=$(awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2; next} FNR==1{for(i=1;i<=NF;i++) if($i=="sample_id") col=i; print; next} {if($col in a) $col=a[$col]; print}' $index_table $table_to_process)
    
    # Execute the AWK script and overwrite the table_to_process
    echo "$awk_script" > $table_to_process
else
    # Scenario 2: The table does not have a header
    
    # Find the matching column number using the first row of the index_table
    first_sample_id=$(awk 'NR==1{print $1}' $index_table)
    matching_col=$(awk -v id="$first_sample_id" '{for(i=1;i<=NF;i++) if($i==id) print i; exit}' $table_to_process)
    
    if [ -z "$matching_col" ]; then
        echo "Could not find matching sample_id in the table to be processed."
        exit 1
    fi
    
    # Generate an AWK script to perform the replacement without assuming headers
    awk_script=$(awk -v col="$matching_col" 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2; next} {if($col in a) $col=a[$col]; print}' $index_table $table_to_process)
    
    # Execute the AWK script and overwrite the table_to_process
    echo "$awk_script" > $table_to_process
fi

