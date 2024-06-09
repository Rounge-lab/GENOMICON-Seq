import gzip
import sys

# Mapping of accession numbers to standard chromosome names
refseq_to_chr = {
    'NC_000001.11': 'chr1', 'NC_000002.12': 'chr2', 'NC_000003.12': 'chr3',
    'NC_000004.12': 'chr4', 'NC_000005.10': 'chr5', 'NC_000006.12': 'chr6',
    'NC_000007.14': 'chr7', 'NC_000008.11': 'chr8', 'NC_000009.12': 'chr9',
    'NC_000010.11': 'chr10', 'NC_000011.10': 'chr11', 'NC_000012.12': 'chr12',
    'NC_000013.11': 'chr13', 'NC_000014.9': 'chr14', 'NC_000015.10': 'chr15',
    'NC_000016.10': 'chr16', 'NC_000017.11': 'chr17', 'NC_000018.10': 'chr18',
    'NC_000019.10': 'chr19', 'NC_000020.11': 'chr20', 'NC_000021.9': 'chr21',
    'NC_000022.11': 'chr22', 'NC_000023.11': 'chrX', 'NC_000024.10': 'chrY'
    # Add other mappings if necessary
}

def process_fasta(input_file, output_file):
    with gzip.open(input_file, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
        write_sequence = False

        for line in infile:
            if line.startswith('>'):
                accession_number = line.split(' ')[0].replace('>', '').strip()
                print(f'Found chromosome: {accession_number}')  

                if accession_number in refseq_to_chr:
                    chr_name = refseq_to_chr[accession_number]
                    print(f'Writing chromosome: {chr_name}') 
                    outfile.write(f'>{chr_name}\n')
                    write_sequence = True
                else:
                    write_sequence = False
            elif write_sequence:
                outfile.write(line)

if __name__ == "__main__":
    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]
    process_fasta(input_file_path, output_file_path)
