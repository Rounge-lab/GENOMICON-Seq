import argparse
import pandas as pd

def filter_column(input_file, column_name, output_file):
    # Load the data from the CSV file, specifying it as a tab-separated file
    # and setting the first column as the index if it doesn't have a header
    data = pd.read_csv(input_file, delimiter='\t', index_col=0)

    # Reset index to turn the index into a regular column (if the index has no name, we'll name it)
    data.reset_index(inplace=True)
    data.columns = ['X'] + list(data.columns[1:])  # Ensuring first column is named 'X' if it wasn't
    
    # Check if the column exists in the DataFrame
    if column_name not in data.columns:
        raise ValueError(f"The specified column '{column_name}' does not exist in the input file.")

    # Extract the specified column along with the 'X' column
    result = data[['X', column_name]]
    
    # Write the resulting DataFrame to a CSV file
    result.to_csv(output_file, index=False)

def main():
    # Setup the argument parser
    parser = argparse.ArgumentParser(description="Extracts a specified column from a tab-separated file and outputs it to a new CSV file.")
    parser.add_argument('--input_sbs', type=str, required=True, help='Path to the input TSV file.')
    parser.add_argument('--column_name', type=str, required=True, help='Name of the column to extract.')
    parser.add_argument('--output_sbs', type=str, required=True, help='Path to the output CSV file.')

    # Parse arguments
    args = parser.parse_args()

    # Run the column filter function
    filter_column(args.input_sbs, args.column_name, args.output_sbs)

if __name__ == "__main__":
    main()

