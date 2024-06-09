from Bio import SeqIO
import csv
import gzip
import argparse

def fasta_to_csv(fasta_file, csv_file):
    if fasta_file.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open

    with open_func(fasta_file, "rt") as fasta, open(csv_file, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["fasta_names", "fasta_lengths"])
        for record in SeqIO.parse(fasta, "fasta"):
            csvwriter.writerow([record.id, len(record.seq)])

def main():
    parser = argparse.ArgumentParser(description='Convert FASTA or gzipped FASTA file to CSV.')
    parser.add_argument('fasta_file', type=str, help='Path to the input FASTA or gzipped FASTA file.')
    parser.add_argument('csv_file', type=str, help='Path to the output CSV file.')

    args = parser.parse_args()
    fasta_to_csv(args.fasta_file, args.csv_file)

if __name__ == "__main__":
    main()
