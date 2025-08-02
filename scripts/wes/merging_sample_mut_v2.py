import pandas as pd
import argparse
import os
import glob

def process_files(mutation_occurrence_file, detailed_mutation_file, output_file):
    try:
        # Check if the occurrence file is empty except for the header
        df1 = pd.read_csv(mutation_occurrence_file)
        if df1.empty:
            print(f"Skipping {mutation_occurrence_file} as it is empty.")
            return  # Skip processing if the file is empty

        # Process the detailed mutation file
        df2 = pd.read_csv(detailed_mutation_file, sep=';')
        
        # Prepare the data from occurrence file
        df1['position'] = df1['mutation'].str.extract('(\d+)').astype(int)
        df1['Variant'] = df1['mutation'].str.extract('([A-Z]$)')
        df1.drop(columns='mutation', inplace=True)
        
        # Prepare the detailed mutation data
        df2['position'] = df2['position'].astype(int)
        df2['Variant'] = df2['Variant'].astype(str)

        # Merge the data
        result = pd.merge(df2, df1, on=['position', 'Variant'], how='left')
        result['occurrence'] = result['occurrence'].fillna(0).astype(int)
        result = result[result['occurrence'] > 0]

        if not result.empty:
            result.to_csv(output_file, index=False)
            print(f"Result saved to {output_file}")
            # Remove input files if output is successfully written
            os.remove(mutation_occurrence_file)
            os.remove(detailed_mutation_file)
            print(f"Input files {mutation_occurrence_file} and {detailed_mutation_file} removed after successful processing.")
        else:
            print(f"No data to output after merging for {mutation_occurrence_file}.")

    except Exception as e:
        print(f"Error processing {mutation_occurrence_file} and {detailed_mutation_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Process mutation files by merging, based on directory input.")
    parser.add_argument("directory", help="Directory containing the CSV files")
    args = parser.parse_args()

    # Find all occurrence files
    occurrence_files = glob.glob(os.path.join(args.directory, '*_mutation_occurrences.csv'))

    for occurrence_file in occurrence_files:
        base_name = os.path.basename(occurrence_file).replace('_mutation_occurrences.csv', '')
        detailed_mutation_file = os.path.join(args.directory, f"{base_name}_mutations_table.csv")
        output_file = os.path.join(args.directory, f"{base_name}_inserted_mutations_overview.csv")

        if os.path.exists(detailed_mutation_file):
            process_files(occurrence_file, detailed_mutation_file, output_file)
        else:
            print(f"No matching detailed mutation table found for {occurrence_file}")

if __name__ == "__main__":
    main()
