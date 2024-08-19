## skimming_pipeline.sh


<Usage: "bash ${BASH_SOURCE[0]} -h [input] [-o output directory] [-t threads]>
Runs nuclear read processing pipeline on a batch of merged and decontaminated reads in reference to a constructed library:


# We will retrieve all the fq sampled by Stein in a folder with the symlinks to run echarvel skimming_pipeline.sh
```bash
module load sratoolkit sra-tools pigz
# Base directory where the CSV file is located
base_dir="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq"
# Path to the CSV file
csv_file="${base_dir}/oedulis_sra.csv"
# Target directory to create symlinks
target_dir="${base_dir}/Stein_raw_fq"
# Create the target directory if it doesn't exist
mkdir -p "$target_dir"
# Navigate to the base directory
cd "$base_dir" || exit
# Read SRA run IDs sampled by Stein Mortensen from the CSV, skip the header line
sra_ids=$(tail -n +2 "$csv_file" | cut -d',' -f1,11 | grep "Stein Mortensen" | cut -d',' -f1)
for sra_id in $sra_ids; do
    # Check if .fastq.gz files already exist for this SRA ID
    if [ -f "${sra_id}_1.fastq.gz" ] && [ -f "${sra_id}_2.fastq.gz" ]; then
        echo "FASTQ files for $sra_id already exist. Creating symlinks."
        
        # Create symlinks in the target directory
        ln -s "${base_dir}/${sra_id}_1.fastq.gz" "${target_dir}/${sra_id}_1.fastq.gz"
        ln -s "${base_dir}/${sra_id}_2.fastq.gz" "${target_dir}/${sra_id}_2.fastq.gz"
    else
        # Download SRA file using prefetch
        prefetch "$sra_id" --output-directory "./"

        # Explicitly define the path of the downloaded SRA file
        sra_file="${sra_id}/${sra_id}.sra"

        # Check if the SRA file was downloaded and exists
        if [ -f "$sra_file" ]; then
            # Convert SRA to FASTQ and compress
            fastq-dump --split-files "$sra_file"
            pigz "${sra_id}"_*.fastq
            rm -f "$sra_file"

            # Create symlinks in the target directory
            ln -s "${base_dir}/${sra_id}_1.fastq.gz" "${target_dir}/${sra_id}_1.fastq.gz"
            ln -s "${base_dir}/${sra_id}_2.fastq.gz" "${target_dir}/${sra_id}_2.fastq.gz"
        else
            echo "Error: SRA file for $sra_id not found."
        fi
    fi
done
echo "All processing and symlink creation complete."
```



# Run on Oedulis Stein raw fastq
```bash
conda activate skimming_echarvel
module purge
module load parallel
bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline_withoutKRANK.sh \
-i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/Stein_raw_fq/ \
-o /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/skmer/echarvel_Stein_fast-skims_results_12.07 -t 18 -p 12 \
-f _1.fastq.gz \
-r _2.fastq.gz
```

# Run on Herring raw fastq
```bash
#Get the raw herring fastq in a new dir
# Define paths
CSV_FILE="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/Clupea_n42_13var.csv"
SOURCE_DIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/rawfastqmodernClupea"
NEW_DIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/fastq"
mkdir -p "$NEW_DIR"
# Extract the cleaned_id values from the CSV file and create symlinks
tail -n +2 "$CSV_FILE" | while IFS=',' read -r cleaned_id _; do
    cleaned_id=$(echo "$cleaned_id" | tr -d '"')  # Remove any surrounding quotes
    for file in "$SOURCE_DIR"/*; do
        filename=$(basename "$file")
        if [[ "$filename" == "$cleaned_id"*".fastq.gz" ]]; then
            ln -s "$file" "$NEW_DIR/$filename"
            echo "Created symlink: $NEW_DIR/$filename -> $file"
        fi
    done
done
echo "Symlinks created successfully."

#skimming_pipeline.sh
conda activate skimming_echarvel
module purge
module load parallel
bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline.sh -i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/fastq -o /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/fastq/echarvel_fast-skims_11.07 -t 18 -p 10
```

# Run on Magpie raw fastq
```bash
conda activate skimming_echarvel
module purge
module load parallel
bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline.sh -i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/done -o /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/echarvel_Stein_fast-skims_results_12.07 -t 18 -p 12
```



# skimming_pipeline_withoutKRANK.sh on herring 
```bash
cp /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline.sh /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline_withoutKRANK.sh

#on herring dataset
conda activate skimming_echarvel
module purge
module load parallel
bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline_withoutKRANK.sh \
-i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/fastq \
-o /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/fastq/echarvel_fast-skims_11.07 -t 18 -p 10 \
-f _1.fastq.gz \
-r _2.fastq.gz
# -> high error-rate ~0.01 and empty dist matrix
```


# skimming_pipeline with only microbe database for KRANK
```bash
#Herring
conda activate skimming_echarvel
module purge
module load parallel
bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline.sh \
-i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/fastq \
-o /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/fastq/echarvel_fast-skims_30.07_onlymicrobe -t 20 -p 10 \
-f _1.fastq.gz \
-r _2.fastq.gz \
-l /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/KRANK/lib_reps_adpt-k29_w35_h13_b16_s8/

#oedulis
bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline.sh \
-i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/Stein_raw_fq/ \
-o /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/skmer/echarvel_Steinoedulis_fast-skims_results_04.08_onlymicrobe -t 18 -p 12 \
-f _1.fastq.gz \
-r _2.fastq.gz \
-l /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/KRANK/lib_reps_adpt-k29_w35_h13_b16_s8/

```

# skimming_pipeline with KRANKv.0.5.0 with only microbe database for KRANK
```bash
#oedulis
conda activate skimming_echarvel
module purge
module load parallel
bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline.sh \
-i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/Stein_raw_fq/ \
-o /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/skmer/echarvel_Steinoedulis_fast-skims_results_06.08_threshold0.05_0.2 -t 18 -p 12 \
-f _1.fastq.gz \
-r _2.fastq.gz 

#herring
bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline.sh \
-i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/fastq \
-o /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/echarvel_fast-skims_07.08_threshold0.05_0.2 -t 20 -p 10 \
-f _1.fastq.gz \
-r _2.fastq.gz 

```