import csv
import sys

def read_primer_names(csv_file):
    with open(csv_file, newline='') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        return {row[0] for row in reader}

def filter_bed_file(bed_file, primer_names, output_file):
    with open(bed_file) as bed, open(output_file, 'w') as out:
        for line in bed:
            parts = line.strip().split('\t')
            if parts[3] in primer_names:
                out.write(line)

def main():
    if len(sys.argv) != 4:
        print("Usage: python filter_bed.py <primers.bed> <primer_set.csv> <output.bed>")
        sys.exit(1)

    bed_file = sys.argv[1]
    csv_file = sys.argv[2]
    output_file = sys.argv[3]

    primer_names = read_primer_names(csv_file)
    filter_bed_file(bed_file, primer_names, output_file)

if __name__ == "__main__":
    main()
