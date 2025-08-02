#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

input_file="$1"
temp_file="$(mktemp)"

awk -F';' '{
    split($1, parts, "_");
    chr = parts[1];
    mg = parts[2];

    if (mg ~ /^NMG/) {
        order = "0";
    } else if (mg ~ /^MG/) {
        split(mg, num, "MG");
        order = num[2];
    } else {
        order = "9999";
    }

    print chr, order, $0
}' "$input_file" | sort -k1,1 -k2,2n | cut -d' ' -f3- > "$temp_file"

mv "$temp_file" "$input_file"
