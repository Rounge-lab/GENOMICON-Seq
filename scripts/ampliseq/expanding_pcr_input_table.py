import pandas as pd
import sys

def process_row(row):
    coords = row['genome_coordinates'].split(',')
    lengths = row['lengths'].split(',')

    if len(coords) > 1:
        new_rows = []
        for i, (coord, length) in enumerate(zip(coords, lengths)):
            new_row = row.copy()
            new_row['genome_coordinates'] = coord
            new_row['lengths'] = length
            new_row['fragment_name'] = f"{row['fragment_name']}-{i+1}"
            new_rows.append(new_row)
        return new_rows
    else:
        return [row]

def main(input_file):
    data = pd.read_csv(input_file, sep=';')

    if not any(',' in str(x) for x in data['genome_coordinates']):
        print("No rows require splitting. File left unchanged.")
        return

    processed_rows = [new_row for _, row in data.iterrows() for new_row in process_row(row)]
    new_data = pd.DataFrame(processed_rows)

    new_data.to_csv(input_file, sep=';', index=False)
    print(f"File '{input_file}' has been successfully processed and overwritten.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>")
        sys.exit(1)

    main(sys.argv[1])
