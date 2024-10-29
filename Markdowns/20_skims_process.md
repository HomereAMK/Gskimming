
# 1.Subsampling raw fastq
```bash
#!/bin/bash

# Variables
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/rawfastqmodernClupea/"
DATAOUTPUT="/projects/mjolnir1/people/sjr729/Skmer_ms/Herring/subsampled_4"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
BBMAP_DIR="/projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/bbmap"

# Create output directory if it doesn't exist
mkdir -p "$DATAOUTPUT"

# Calculate genome size
genome_size=$(grep -v "^>" "$GENOME" | awk '{total += length($0)} END {print total}')

# Desired coverage
coverage=4

# Average read length (adjust if your reads are different)
avg_read_length=150

# Calculate desired total bases
desired_bases=$((genome_size * coverage))

# Number of samples (assuming paired-end reads)
num_samples=$(ls "$DATAINPUT"/*_1.fastq.gz | wc -l)

# Desired reads per sample
reads_per_sample=$((desired_bases / avg_read_length / num_samples))

echo "Genome size: $genome_size bp"
echo "Desired total bases: $desired_bases bp"
echo "Number of samples: $num_samples"
echo "Reads per sample per sample: $reads_per_sample"

# Subsample FASTQ files
for file1 in "$DATAINPUT"/*_1.fastq.gz; do
    # Determine the base name
    base=$(basename "$file1" _1.fastq.gz)
    file2="${DATAINPUT}/${base}_2.fastq.gz"

    # Subsampled output files in DATAOUTPUT
    subsampled_file1="${DATAOUTPUT}/subsampled4x_${base}_1.fastq.gz"
    subsampled_file2="${DATAOUTPUT}/subsampled4x_${base}_2.fastq.gz"

    # Check if both subsampled files already exist
    if [[ -f "$subsampled_file1" && -f "$subsampled_file2" ]]; then
        echo "Subsampled files already exist for $base. Skipping subsampling."
        continue
    fi

    # Check if the reverse read file exists
    if [[ -f "$file2" ]]; then
        echo "Subsampling $base"
        "$BBMAP_DIR"/reformat.sh \
            in1="$file1" in2="$file2" \
            out1="$subsampled_file1" out2="$subsampled_file2" \
            samplereadstarget="$reads_per_sample" \
            overwrite=false
    else
        echo "Reverse read file not found for $file1. Skipping."
    fi
done



```

# 2.1  Test bbmap_pipeline_oct24.sh
```bash
INPUTDIR=/projects/mjolnir1/people/sjr729/Skmer_ms/Herring/subsampled_4
cd $INPUTDIR
/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/bbmap_pipeline_oct24.sh \
subsampled4x_ERR9709203_1.fastq.gz \
subsampled4x_ERR9709203_2.fastq.gz \
subsampled4x_ERR9709203_merged.fq.gz
```





# 2.2 New preprocess without bbmap merge step
```bash
#For Herring
conda activate skimming_echarvel
module load parallel
OUTPUTDIR=/projects/mjolnir1/people/sjr729/Skmer_ms/Herring
DATAINPUT="/projects/mjolnir1/people/sjr729/Skmer_ms/Herring/subsampled_4"

cd $OUTPUTDIR

/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/skims_processing_pipeline_oct24.sh -x $DATAINPUT -r 30 -f 30 > $OUTPUTDIR/skims_14oct24.log 2>&1


#For Oedulis, no subsampling
conda activate skimming_echarvel
module load parallel
OUTPUTDIR=/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/fastqoct24"

cd $OUTPUTDIR

/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/skims_processing_pipeline_oct24.sh -x $DATAINPUT -r 10 -f 10 > $OUTPUTDIR/skims_oedulis_25oct24_3.log


```