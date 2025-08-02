import pandas as pd
import os
import sys
import glob

def safe_load_csv(file_path):
    """Safely loads a CSV file, ensuring it is not empty before reading."""
    if os.stat(file_path).st_size == 0:
        print(f"Warning: The file {file_path} is empty.")
        return pd.DataFrame()
    else:
        return pd.read_csv(file_path)

def update_overview_table(directory_path):
    # Debug: Print the input directory path
    print(f"Processing directory: {directory_path}")

    # Construct the pattern for the main data table and the log file
    main_table_pattern = os.path.join(directory_path, "*-real_inserted_mutations.csv")
    log_file_pattern = os.path.join(directory_path, "*_no_introduced_true_mutations.log")

    # Find files matching the patterns
    main_table_files = glob.glob(main_table_pattern)
    log_file_exists = glob.glob(log_file_pattern)


    # Check for log files and absence of data files
    if log_file_exists or not main_table_files:
        print("No mutations were inserted, skipping processing.")
        return

    main_table_path = main_table_files[0]

    
    main_df = safe_load_csv(main_table_path)
    if main_df.empty:
        print("No data found in the main file. No processing will be performed.")
        os.remove(main_table_path)
        print(f"Main input file {main_table_path} has been removed as it was empty.")
        return

    main_df.columns = ['chromosome', 'mutation', 'value']

    unique_identifiers = main_df['chromosome'].unique()
    overview_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(main_table_path)))))

    for identifier in unique_identifiers:
        overview_table_name = f"{identifier}_inserted_mutations_overview.csv"
        overview_table_path = os.path.join(overview_path, overview_table_name)
        print(f"Looking for overview table at: {overview_table_path}")  # Debug

        if os.path.exists(overview_table_path):
            overview_df = pd.read_csv(overview_table_path)
        else:
            print(f"Overview table not found for {identifier}, skipping.")
            continue

        matched_rows = pd.DataFrame(columns=overview_df.columns)

        for _, row in main_df.iterrows():
            numeric_part = ''.join(filter(str.isdigit, row['mutation']))
            variant_letter = ''.join(filter(str.isalpha, row['mutation']))

            mask = (overview_df['position'] == int(numeric_part)) & (overview_df['Variant'] == variant_letter)
            overview_df.loc[mask, 'occurrence'] = row['value']
            matched_rows = matched_rows._append(overview_df[mask])

        output_file_path = os.path.join(os.path.dirname(main_table_path), f"{identifier}_seq_mutations_overview.csv")
        matched_rows.to_csv(output_file_path, index=False)
        print(f"Updated table saved to {output_file_path}")

    os.remove(main_table_path)
    print(f"Main input file {main_table_path} has been removed successfully.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <directory_with_main_table>")
        sys.exit(1)
    directory_with_main_table = sys.argv[1]
    update_overview_table(directory_with_main_table)
