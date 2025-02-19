
```bash
# Variables
DATAINPUT="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/HC_fastq/fastq"
DATAOUTPUT="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/HC_fastq/subsampled_4"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
BBMAP_DIR="/projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/bbmap"

cd $DATAINPUT
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
num_samples=$(ls "$DATAINPUT"/*_1.fq.gz | wc -l)

# Desired reads per sample
reads_per_sample=$((desired_bases / avg_read_length ))

echo "Genome size: $genome_size bp"
echo "Desired total bases: $desired_bases bp"
echo "Number of samples: $num_samples"
echo "Reads per sample per sample: $reads_per_sample"

# Subsample FASTQ files
for file1 in "$DATAINPUT"/*_1.fq.gz; do
    # Determine the base name
    base=$(basename "$file1" _1.fq.gz)
    file2="${DATAINPUT}/${base}_2.fq.gz"

    # Subsampled output files in DATAOUTPUT
    subsampled_file1="${DATAOUTPUT}/subsampled4x_${base}_1.fq.gz"
    subsampled_file2="${DATAOUTPUT}/subsampled4x_${base}_2.fq.gz"

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

#
```bash
conda activate skimming_echarvel
module load parallel
OUTPUTDIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/HC_fastq/subsampled_4/skimpreprocpip_new"
DATAINPUT="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/HC_fastq/subsampled_4/fastq"

mkdir -p $OUTPUTDIR
cd $OUTPUTDIR
/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/skims_processing_pipeline_oct24.sh -x $DATAINPUT -r 30 -f 30 > $OUTPUTDIR/HC_oedulis_skims_11dec24.log 2>&1
```



# repeat spectra GenomeFR
```bash
OUT_DIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/HC_fastq/subsampled_4/Skmer2"

genome_FR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/Tanguy/fileOegenome10scaffoldC3G.fasta"
GENOMEUK="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/ncbi_dataset/data/GCF_947568905.1/GCF_947568905.1_xbOstEdul1.1_genomic.fna"

conda activate skimming_echarvel
module load parallel
threads=30

cd $OUT_DIR
respect -i ${GENOMEUK} --threads ${threads}

#paste the estimated genome length to estimated-spectra.txt as a last column
#paste <(cat estimated-spectra.txt) <(cut -f5 estimated-parameters.txt) > respect-FR-reference.txt

```


# Skmer1
```bash
INPUT_DIR=/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/HC_fastq/subsampled_4/skimpreprocpip/skims_processing_pipeline_oct24/kraken/
OUT_DIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/HC_fastq/subsampled_4/Skmer1"
conda activate skimming_echarvel
module load parallel

cd $OUT_DIR
/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/run_skmer.sh -i $INPUT_DIR -t 30 -p 10 -o $OUT_DIR/skmer1_library
```


# Repeat spectra with respect on skmer1 library
```bash
NUM_THREADS=10
OUTPUT_DIRECTORY="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/HC_fastq/subsampled_4/"

cd $OUTPUT_DIRECTORY

# Create separate directories for .hist files and hist_info.txt
mkdir --parents "${OUTPUT_DIRECTORY}/respect/hist_files/"
echo -e "Input\tread_length" > "${OUTPUT_DIRECTORY}/respect/hist_info.txt"

# Adjusted the path to match your directory structure
for directory in $(find "${OUTPUT_DIRECTORY}/Skmer1/skmer1_library/skmer_library" -maxdepth 1 -mindepth 1 -type d); do
    file=${directory##*/}
    
    # Ensure get_read_length function is available or replace with appropriate command
    read_len=$(grep 'read_length' "${directory}/${file}.dat" | awk '{print $2}')
    
    echo -e "${file}.hist\t${read_len}" >> "${OUTPUT_DIRECTORY}/respect/hist_info.txt"
    ln --symbolic "$(realpath "${directory}/${file}.hist")" "${OUTPUT_DIRECTORY}/respect/hist_files/"
done

# Now run respect, pointing to the hist_files directory
respect --input-directories "${OUTPUT_DIRECTORY}/respect/hist_files/" \
        --info-file "${OUTPUT_DIRECTORY}/respect/hist_info.txt" \
        --output-directory "${OUTPUT_DIRECTORY}/respect/output/" \
        --spectra-output-size 50 \
        --iterations 1000 \
        --threads "${NUM_THREADS}"


cd $OUTPUT_DIRECTORY
#paste the estimated genome length to estimated-spectra.txt as a last column
paste <(cat respect/output/estimated-spectra_1.txt) <(cut -f5 respect/output/estimated-parameters_1.txt) > respect/output/respect-hist-reference.txt
```



# Skmer2 UK and FR genome
```bash
/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/HC_fastq/subsampled_4/skimpreprocpip/skims_processing_pipeline_oct24/kraken/
conda activate skimming_echarvel



module load parallel
RESPECT_SPECTRA_UK=/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/respect-reference.txt
INPUT_DIR=/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/HC_fastq/subsampled_4/skimpreprocpip/skims_processing_pipeline_oct24/kraken/
OUT_DIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/HC_fastq/subsampled_4/Skmer2"
GENOMEUK="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/ncbi_dataset/data/GCF_947568905.1/GCF_947568905.1_xbOstEdul1.1_genomic.fna"
genome_FR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/Tanguy/fileOegenome10scaffoldC3G.fasta"
Respect_spectra_hist="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/HC_fastq/subsampled_4/respect/output/respect-hist-reference.txt"

cd $OUT_DIR
#python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/skmer_new-err.py --debug reference $INPUT_DIR -r $RESPECT_SPECTRA_UK -l $OUT_DIR/skmer2_library_UKgenome_respect_ref/ -p 1

#skmer distance $OUT_DIR/skmer2__FRgenome_ref_library_v3 -t 2 -o FR_genome_ref_library-dist-mat

#python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/skmer_new-err.py --debug reference $INPUT_DIR -r $Respect_spectra_hist  -l $OUT_DIR/skmer2_library_respect_spectra_hist_ref/ -p 1

python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/skmer_new-err.py --debug reference $INPUT_DIR -r $GENOMEUK -l $OUT_DIR/skmer2__UKgenome_ref_library_v3/ -p 1

python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/skmer_new-err.py --debug reference $INPUT_DIR -r $genome_FR -l $OUT_DIR/skmer2__FRgenome_ref_library_v3/ -p 1



# Skimpreprocess raw HCFastq

conda activate skimming_echarvel
module load parallel
OUTPUTDIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/HC_fastq/Raw_skimpreprocpip"
DATAINPUT="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/HC_fastq/fastq"

mkdir -p $OUTPUTDIR
cd $OUTPUTDIR

bash /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/skims_processing_pipeline_oct24.sh -x $DATAINPUT -r 30 -f 30