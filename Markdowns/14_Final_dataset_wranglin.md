## Oedulis
## Sampled skim_processed .fq in fastq/
```bash
#!/bin/bash

# Define the base directory and CSV file path
base_dir="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/"
csv_file="${base_dir}/oedulis_sra.csv"
destination_dir="${base_dir}/oedulis_RawFastq_14.08.24"

# Create the destination directory if it does not exist
mkdir -p "$destination_dir"

# List of identifiers to grep for
identifiers=("PONT" "TOLL" "TRAL" "RYAN" "NISS" "HYPP" "BUNN" "AGAB" "OSTR" "VAGS")

# Loop over each identifier
for identifier in "${identifiers[@]}"; do
    echo "Processing identifier: $identifier"
    
    # Get SRA run IDs for the current identifier
    sra_ids=$(tail -n +2 "$csv_file" | grep "$identifier" | cut -d',' -f1)
    
    # Loop over each SRA ID
    for sra_id in $sra_ids; do
        # Check if the matching fastq.gz files already exist in the destination directory
        if [ -f "${destination_dir}/${sra_id}_1.fastq.gz" ] && [ -f "${destination_dir}/${sra_id}_2.fastq.gz" ]; then
            echo "Skipped: Files already exist for $sra_id"
        else
            # Check if the SRA files exist in the base directory
            if [ -f "${sra_id}_1.fastq.gz" ] && [ -f "${sra_id}_2.fastq.gz" ]; then
                # Copy the matching fastq.gz files to the destination directory
                cp "${sra_id}_1.fastq.gz" "${destination_dir}/"
                cp "${sra_id}_2.fastq.gz" "${destination_dir}/"
                echo "Copied: ${sra_id}_1.fastq.gz and ${sra_id}_2.fastq.gz to $destination_dir"
            else
                echo "File not found: ${sra_id}_1.fastq.gz or ${sra_id}_2.fastq.gz"
            fi
        fi
    done
done
```