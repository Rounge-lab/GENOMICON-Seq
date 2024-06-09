import pandas as pd
import sys
from itertools import product

# Define a dictionary to map ambiguous bases to normal bases
ambiguous_bases = {
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

def replace_ambiguous_bases(seq):
    ambiguous_bases_in_seq = [(i, base) for i, base in enumerate(seq) if base in ambiguous_bases]

    combinations = product(*(ambiguous_bases[base] for _, base in ambiguous_bases_in_seq))

    for combination in combinations:
        new_seq = list(seq)  
        for (i, _), replacement in zip(ambiguous_bases_in_seq, combination):
            new_seq[i] = replacement
        yield ''.join(new_seq)  

def convert_csv_to_fasta(csv_file, fasta_file):
    df = pd.read_csv(csv_file)
    with open(fasta_file, 'w') as f:
        for index, row in df.iterrows():
            for new_seq in replace_ambiguous_bases(row['seq']):
                f.write('>' + row['name'] + '\n')
                f.write(new_seq + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_csv_to_fasta.py <input_csv_file> <output_fasta_file>")
        sys.exit(1)

    csv_file = sys.argv[1]
    fasta_file = sys.argv[2]
    convert_csv_to_fasta(csv_file, fasta_file)