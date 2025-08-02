import pandas as pd
import argparse
import os
import glob
from concurrent.futures import ThreadPoolExecutor

def process_file(file_path, output_directory):
    try:
        # Prepare output file path for occurrences
        base_name = os.path.basename(file_path)
        output_file = os.path.join(output_directory, base_name.replace('sample_table', 'mutation_occurrences'))
        count_output_file = os.path.join(output_directory, base_name.replace('sample_table', 'mutation_counts'))

        # Read and process mutations
        df = pd.read_csv(file_path, sep=',', usecols=['genome_name', 'variant_position', 'copy_number'])
        df['copy_number'] = pd.to_numeric(df['copy_number'], errors='coerce')
        df.dropna(subset=['copy_number'], inplace=True)

        # Create a DataFrame for mutation counts per genome
        counts_df = df.copy()
        counts_df['mutation_number'] = counts_df['variant_position'].apply(lambda x: 0 if x == 'No mutation' else len(x.split(',')))
        counts_df = counts_df[['genome_name', 'mutation_number', 'copy_number']]
        counts_df.to_csv(count_output_file, index=False)
        print(f"Mutation count data saved to '{count_output_file}'")

        # Process mutation occurrences
        df = df[df['variant_position'] != 'No mutation']
        mutations_expanded = df['variant_position'].str.split(',', expand=True).stack().reset_index(level=1, drop=True)
        mutations_expanded.name = 'mutation'
        df = df.drop(columns='variant_position').join(mutations_expanded)
        df['copy_number'] = df['copy_number'].astype(int)
        result_df = df.groupby('mutation')['copy_number'].sum().reset_index(name='occurrence')
        result_df.to_csv(output_file, index=False)

        # Remove the original file upon successful processing
        os.remove(file_path)
        print(f"Processing and removal complete. Output saved to '{output_file}' and original file removed.")
    except Exception as e:
        print(f"Failed to process {file_path} due to {e}.")

def main():
    parser = argparse.ArgumentParser(description="Process multiple genome mutation data files in parallel, generate mutation counts, and remove originals.")
    parser.add_argument("directory", help="Directory containing the input CSV files")
    parser.add_argument("--cores", type=int, default=4, help="Number of cores to use for parallel processing")
    
    args = parser.parse_args()
    input_files = glob.glob(os.path.join(args.directory, '*_sample_table.csv'))
    
    # Adjust the number of cores if more than the number of files
    num_cores = min(args.cores, len(input_files))
    
    with ThreadPoolExecutor(max_workers=num_cores) as executor:
        # Submit processing tasks
        for file_path in input_files:
            executor.submit(process_file, file_path, args.directory)

if __name__ == "__main__":
    main()
