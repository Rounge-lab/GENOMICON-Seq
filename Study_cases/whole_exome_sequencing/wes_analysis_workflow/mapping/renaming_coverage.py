import pandas as pd
import os
import sys

def update_sample_ids(index_table_path, table_to_process_path):
    # Load index table into a dictionary
    index_df = pd.read_csv(index_table_path, sep='\t', dtype=str)
    sample_id_to_index = dict(zip(index_df['sample_id'], index_df['sample_index']))

    # Define a function to process each chunk
    def process_chunk(chunk):
        chunk['sample_id'] = chunk['sample_id'].apply(lambda x: sample_id_to_index.get(x, x))
        return chunk

    # Temporary file path
    temp_output_path = table_to_process_path + '.tmp'

    # Process in chunks
    chunk_size = 10**6  # Adjust based on your system's capabilities
    reader = pd.read_csv(table_to_process_path, sep='\t', dtype=str, chunksize=chunk_size)
    first_chunk = True
    for chunk in reader:
        processed_chunk = process_chunk(chunk)
        if first_chunk:
            processed_chunk.to_csv(temp_output_path, mode='w', index=False, sep='\t')
            first_chunk = False
        else:
            processed_chunk.to_csv(temp_output_path, mode='a', index=False, sep='\t', header=False)

    # Replace the original file with the processed file
    os.replace(temp_output_path, table_to_process_path)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: script.py <index_table_path> <table_to_process_path>")
        sys.exit(1)
    update_sample_ids(sys.argv[1], sys.argv[2])
