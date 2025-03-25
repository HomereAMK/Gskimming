# 1. sftp files


#!/bin/bash
REMOTE_DIR="/scratch03/seenisftp/upload_plant/"
LOCAL_DIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Ragweed/bam"
SFTP_USER="seenisftp"
REMOTE_HOST="calab-fe.ucsd.edu"

# Create a list of remote files
sftp $SFTP_USER@$REMOTE_HOST <<EOF | grep ".bam$" > remote_files.txt


# 2. Genome formatting

# Load necessary modules
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18 jdk/1.8.0_291 picard/2.27.5 parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15 python/3.11.3 openjdk/17.0.8 gatk/3.8 bedtools bbmap mosdepth
# Var
GENOME=/projects/mjolnir1/people/sjr729/Skmer_ms/Ragweed/genome/ncbi_dataset/data/GCA_032199755.1/GCA_032199755.1_A.artemisiifolia_dipasm_2023_h1_p_genomic.fna
# Create dependencies for the genome
samtools faidx $GENOME
bwa index $GENOME
picard CreateSequenceDictionary R=$GENOME O=${GENOME%.fna}.dict

# 3. Create table for loco-pipe

```bash
bam_dir="/projects/mjolnir1/people/sjr729/Skmer_ms/Ragweed/bam"
output_file="/projects/mjolnir1/people/sjr729/Skmer_ms/Ragweed/sample_table_ragweed.tsv"
csv_file="/projects/mjolnir1/people/sjr729/Skmer_ms/Ragweed/Master_list_Ragweed_modern.csv"
species="Ambrosia_artemisiifolia"

echo -e "sample_name\tspecies\tpopulation\tbam" > "$output_file"

for bam_file in "$bam_dir"/*.bam; do
    filename=$(basename "$bam_file")
    sample_name="${filename%%.*}"  # Extract text before the first '.'

    # Use Python to find the population for the sample_name
    population=$(python3 -c "
import csv
with open('$csv_file', newline='', encoding='utf-8') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if row['Individual ID'].strip('\"') == '$sample_name':
            print(row['Spatial group assignment'].strip('\"'))
            break
    else:
        print('Unknown')
")

    if [ -z "$population" ]; then
        population="Unknown"
    fi

    # Remove any spaces, commas, or tabs from the population variable
    population_clean=$(echo "$population" | tr -d ' ,\t')

    echo -e "${sample_name}\t${species}\t${population_clean}\t${bam_file}" >> "$output_file"
done
```


# make the chr file 
```
# Using cut:
cut -f 1 genome/ragweed-dipasm-hap1.fasta.fai > loco-pipe-git_1e-3/docs/chr_ragweed-dipasm-hap1.tsv
awk '$1 ~ /^h1s([1-9]|1[0-8])$/ { print }' loco-pipe-git_1e-3/docs/chr_ragweed-dipasm-hap1.tsv > loco-pipe-git_1e-3/docs/chr_ragweed-dipasm-hap1_subset_h1s1-18.fai

```


# Set up the file structure.
```bash

SOFTWARE_DIR=/projects/mjolnir1/people/sjr729/
BASEDIR=/projects/mjolnir1/people/sjr729/Skmer_ms/Ragweed
mkdir -p $BASEDIR
cd $BASEDIR
mkdir docs
mkdir config
conda activate loco-pipe-git

#/projects/mjolnir1/people/sjr729/Skmer_ms/Ragweed
#/genome/ncbi_dataset/data//GCA_032199755.1/GCA_032199755.1_A.artemisiifolia_dipasm_2023_h1_p_genomic.fna
#/projects/mjolnir1/people/sjr729/Skmer_ms/Ragweed/docs/sample_table_ragweed.tsv
#/projects/mjolnir1/people/sjr729/Skmer_ms/Ragweed/docs/chr_list.tsv

conda activate loco-pipe-git
module purge
SOFTWARE_DIR=/projects/mjolnir1/people/sjr729/
BASEDIR=/projects/mjolnir1/people/sjr729/Skmer_ms/Ragweed/loco-pipe-git_1e-3
cd $SOFTWARE_DIR
snakemake \
--use-conda \
--conda-frontend mamba \
--directory $BASEDIR \
--rerun-triggers mtime \
--scheduler greedy \
--printshellcmds \
--snakefile $SOFTWARE_DIR/loco-pipe/workflow/pipelines/loco-pipe.smk \
--cores 10 --rerun-incomplete --keep-going 


