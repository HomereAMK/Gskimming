import csv
import sys

def process_csv_to_tsv(input_csv_file, output_tsv_file):
    with open(input_csv_file, mode='r') as infile:
        # Adjust the delimiter here if it's not a comma
        reader = csv.reader(infile, delimiter='\t')

        data_dict = {}

        for row in reader:
            # Skip rows that don't have exactly 3 columns
            if len(row) != 3:
                print(f"Skipping row: {row}")
                continue

            key, variable, value = row

            if key not in data_dict:
                data_dict[key] = {}
            data_dict[key][variable] = value

    with open(output_tsv_file, mode='w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')

        for key, values in data_dict.items():
            for variable, value in values.items():
                writer.writerow([key, variable, value])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python csv_to_tsv_stats.py <input_csv_file> <output_tsv_file>")
        sys.exit(1)

    input_csv_file = sys.argv[1]
    output_tsv_file = sys.argv[2]

    process_csv_to_tsv(input_csv_file, output_tsv_file)

(tutorial) [sjr729@mjolnirhead01fl testClupea]$ cat badspecimens_identifier.py

import sys
import os

def move_bad_specimens(tsv_file, source_dir, destination_dir):
    os.makedirs(destination_dir, exist_ok=True)
    specimen_data = {}

    with open(tsv_file, 'r') as file:
        for line in file:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue

            file_base, criterion, value = parts
            filename = os.path.join(source_dir, file_base.strip())

            if filename not in specimen_data:
                specimen_data[filename] = {'error_rate': None, 'read_length': None}
            if 'error_rate' in criterion:
                specimen_data[filename]['error_rate'] = value
            elif 'read_length' in criterion:
                specimen_data[filename]['read_length'] = value

    moved_files_count = 0
    for filename, data in specimen_data.items():
        error_rate = data['error_rate']
        read_length = data['read_length']

        try:
            if (error_rate == 'NA' or float(error_rate) > 0.02) or (read_length != 'NA' and int(read_length) < 70):
                if os.path.exists(filename):
                    dest_file = os.path.join(destination_dir, os.path.basename(filename))
                    print(f"Moving {filename} to {dest_file}")
                    os.rename(filename, dest_file)
                    moved_files_count += 1
                else:
                    print(f"File not found: {filename}")
        except Exception as e:
            print(f"Error moving file {filename}: {e}")

    print(f"{moved_files_count} files moved to {destination_dir}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python badspecimens_identifier.py input.tsv source_directory destination_directory")
        sys.exit(1)

    input_tsv = sys.argv[1]
    source_dir = sys.argv[2]
    destination_dir = sys.argv[3]

    move_bad_specimens(input_tsv, source_dir, destination_dir)
