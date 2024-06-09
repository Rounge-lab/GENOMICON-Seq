#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <csv1> <csv2> <csv3> <output_csv>"
    exit 1
fi

csv1=$1
csv2=$2
csv3=$3
output_csv=$4

awk -F';' '{
    # Find the position of the fourth underscore
    split($1, parts, "_");
    prefix = parts[1];
    for (i = 2; i <= 4; i++) {
        prefix = prefix "_" parts[i];
    }
    suffix = substr($1, index($1, parts[5]));
    print prefix ";" suffix ";" $2
}' $csv1 > temp1.csv

awk -F';' 'BEGIN{OFS=FS} {gsub(" ","_",$5)} 1' $csv2 | awk -F';' 'NR==FNR{refs[$1";"$2]; next} ($1";"$5 in refs)' temp1.csv - | awk -F';' '{print $1 ";" $5}' > temp2.csv

awk -F';' 'NR==FNR{info[$1]=$2; next} $1 in info {print $0 ";" info[$1]}' $csv3 temp2.csv > "$output_csv"

rm temp1.csv temp2.csv

echo "Processing complete. Output saved to $output_csv"

