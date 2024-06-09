#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_filename> <output_filename>"
    exit 1
fi

INPUT_FILE=$1
OUTPUT_FILE=$2

awk -F';' 'NR > 1 { 
    split($4, coords, ":"); 
    diff = coords[2] - coords[1] + 1; 
    print $1 "_" $2 ";" diff 
}' $INPUT_FILE > $OUTPUT_FILE

