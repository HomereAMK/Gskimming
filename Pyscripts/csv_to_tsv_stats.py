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
