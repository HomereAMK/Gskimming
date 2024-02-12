import csv
import subprocess
import sys

# Check if a command-line argument was provided
if len(sys.argv) < 2:
    print("Usage: python script.py <path_to_your_csv_file>")
    sys.exit(1)

# The first command-line argument after the script name is the CSV file path
csv_file = sys.argv[1]

with open(csv_file, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        sample_id = row[0]  # Assuming the sample ID is in the first column
        samplerate = row[5]  # Assuming the samplerate is in the sixth column

        # Ensure samplerate is a float and handle any potential errors
        try:
            samplerate = float(samplerate)
        except ValueError:
            # Handle the error or skip this row
            continue

        # Construct the reformat.sh command
        command = f'reformat.sh in1=done/{sample_id}_1.fastq.gz in2=done/{sample_id}_2.fastq.gz out1=./{sample_id}_4x_1.fastq.gz out2=./{sample_id}_4x_2.fastq.gz samplerate={samplerate:.6f}'

        # Execute the command
        try:
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while processing {sample_id}: {e}")
