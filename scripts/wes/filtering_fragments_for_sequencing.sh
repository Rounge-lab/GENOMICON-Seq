#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <first_file.csv> <second_file.csv> <output_file.csv>"
    exit 1
fi

first_file="$1"
second_file="$2"
output_file="$3"

awk -F ';' '{ gsub("_No_mutations", ""); print $1 }' "$first_file" > names_temp.txt

awk -F ';' 'NR==FNR { names[$1]; next } $1 in names' names_temp.txt "$second_file" > "$output_file"

rm names_temp.txt
rm "$second_file"

echo "Filtering complete. Output is in $output_file"
