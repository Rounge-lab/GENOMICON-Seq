import pandas as pd
import os
import sys

def safe_load_csv(file_path):
    """Safely loads a CSV file, ensuring it is not empty before reading."""
    if os.stat(file_path).st_size == 0:  # Check if the file size is 0 bytes
        print(f"Warning: The file {file_path} is empty.")
        return pd.DataFrame()  # Return an empty DataFrame if file is empty
    else:
        return pd.read_csv(file_path, header=None)

def update_overview_table(main_table_path):
    # Load the main table without headers, safely
    main_df = safe_load_csv(main_table_path)
    if main_df.empty:
        print("No data found in the main file. No processing will be performed.")
        os.remove(main_table_path)  # Remove the empty file
        print(f"Main input file {main_table_path} has been removed as it was empty.")
        return

    # Set column names manually since the file does not have headers
    main_df.columns = ['chromosome', 'mutation', 'value']

    # Extract unique identifiers from the first column
    unique_identifiers = main_df['chromosome'].unique()

    # Deduce the directory path for the overview tables
    overview_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(main_table_path))))

    # Process each unique identifier
    for identifier in unique_identifiers:
        # Load the corresponding overview table
        overview_table_name = f"{identifier}_inserted_mutations_overview.csv"
        overview_table_path = os.path.join(overview_path, overview_table_name)
        overview_df = pd.read_csv(overview_table_path)

        # Create an empty DataFrame to store matched rows
        matched_rows = pd.DataFrame(columns=overview_df.columns)

        # Iterate over the main_df to update the overview_df
        for _, row in main_df.iterrows():
            # Split numeric part and variant letter from 'mutation'
            numeric_part = ''.join(filter(str.isdigit, row['mutation']))
            variant_letter = ''.join(filter(str.isalpha, row['mutation']))

            # Find matching rows and update 'occurrence'
            mask = (overview_df['position'] == int(numeric_part)) & (overview_df['Variant'] == variant_letter)
            overview_df.loc[mask, 'occurrence'] = row['value']
            matched_rows = matched_rows._append(overview_df[mask])

        # Define output file path
        output_file_path = os.path.join(os.path.dirname(main_table_path), f"{identifier}_seq_mutations_overview.csv")

        # Save only the matched rows to the output file
        matched_rows.to_csv(output_file_path, index=False)
        print(f"Updated table saved to {output_file_path}")

    # Remove the main input file after processing
    os.remove(main_table_path)
    print(f"Main input file {main_table_path} has been removed successfully.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_main_table>")
        sys.exit(1)
    path_to_main_table = sys.argv[1]
    update_overview_table(path_to_main_table)
