# 1. 
```bash
NUM_THREADS=10
OUTPUT_DIRECTORY=…

cd $OUTPUT_DIRECTORY

# Create separate directories for .hist files and hist_info.txt
mkdir --parents "${OUTPUT_DIRECTORY}/respect/hist_files/"
echo -e "Input\tread_length" > "${OUTPUT_DIRECTORY}/respect/hist_info.txt"

# Adjusted the path to match your directory structure
for directory in $(find "${OUTPUT_DIRECTORY}/skmer1/skmer1_library/skmer_library" -maxdepth 1 -mindepth 1 -type d); do
    file=${directory##*/}
    
    # Ensure get_read_length function is available or replace with appropriate command
    read_len=$(grep 'read_length' "${directory}/${file}.dat" | awk '{print $2}')
    
    echo -e "${file}.hist\t${read_len}" >> "${OUTPUT_DIRECTORY}/respect/hist_info.txt"
    ln --symbolic "$(realpath "${directory}/${file}.hist")" "${OUTPUT_DIRECTORY}/respect/hist_files/"
done

# Now run respect, pointing to the hist_files directory
respect --input-directories "${OUTPUT_DIRECTORY}/respect/hist_files/" \
        --info-file "${OUTPUT_DIRECTORY}/respect/hist_info.txt" \
        --output-directory "${OUTPUT_DIRECTORY}/" \
        --spectra-output-size 50 \
        --iterations 1000 \
        --threads "${NUM_THREADS}"



cd /projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis
#paste the estimated genome length to estimated-spectra.txt as a last column
paste <(cat estimated-spectra.txt) <(cut -f5 estimated-parameters.txt) > respect-reference.txt


conda activate skimming_echarvel
module load parallel
RESPECT_SPECTRA=/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/respect-reference.txt
INPUT_DIR=/projects/mjolnir1/people/sjr729//tutorial/skimming_scripts/oedulis/fastq/skims_processing_pipeline/kraken/

GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/ncbi_dataset/data/GCF_947568905.1/GCF_947568905.1_xbOstEdul1.1_genomic.fna"

python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/skmer_new-err.py --debug reference $INPUT_DIR -r $RESPECT_SPECTRA -l ./skmer2_library_respect_ref/ -p 2

python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/skmer_new-err.py --debug reference $INPUT_DIR -r $GENOME -l ./skmer2_library_v3/ -p 2
```


# 2. foldering
```bash
# Step 1: Set Variables
OUTPUT_DIR=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/skims_processing_pipeline/foldering
GENOME=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/ncbi_dataset/data/GCF_947568905.1/GCF_947568905.1_xbOstEdul1.1_genomic.fna
PROCESS_SCRIPT=process_dir.sh

# Step 2: Move Each .fq File into Its Own Directory Named dir1 to dir268
#cd "$OUTPUT_DIR"
#N=1
#for DIR in dir{1..14}; do
#    if [ -d "$DIR" ]; then
#        for FILE in "$DIR"/*.fq; do
#            NEW_DIR="$OUTPUT_DIR/dir$N"
#            mkdir -p "$NEW_DIR"
#            mv "$FILE" "$NEW_DIR/"
#            echo "Moved $FILE to $NEW_DIR/"
#            N=$((N + 1))
#        done
#        rmdir "$DIR"
#    else
#        echo "Directory $DIR does not exist."
#    fi
#done
#echo "Total directories created: $((N - 1))"

# Step 3: Create the Processing Script
cd "$OUTPUT_DIR"
echo '#!/bin/bash
DIR=$1
python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/skmer_new-err.py \
--debug reference '"$OUTPUT_DIR"'/$DIR \
-r '"$GENOME"' \
-l '"$OUTPUT_DIR"'/skmer2_library_$DIR/ \
-p 2' > "$PROCESS_SCRIPT"
chmod +x "$PROCESS_SCRIPT"

# Step 4: Load Modules and Activate Environment
module load parallel
conda activate skimming_echarvel

# Step 5: Run the Processing Script Using GNU Parallel
DIRS=(dir*)
parallel -j 10 ./process_dir.sh ::: "${DIRS[@]}"


# Step 6: Combine the Library Directories
mkdir -p skmer2_library_combined
for DIR in dir*; do
    LIB_DIR="skmer2_library_$DIR"
    if [ -d "$LIB_DIR" ]; then
        cp -r "$LIB_DIR/"* skmer2_library_combined/
        echo "Copied contents of $LIB_DIR to skmer2_library_combined/"
    else
        echo "Library directory $LIB_DIR does not exist."
    fi
done

# Optional: Verify Combined Library Contents
echo "Total files in combined library: $(ls skmer2_library_combined/ | wc -l)"



python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/skmer_new-err.py --debug reference $INPUT_DIR -r $GENOME -l ./skmer2_library_v3/ -p 2

```


# Using skmer1 library for skmer2
```bash

module load parallel
conda activate skimming_echarvel


INPUT_DIR=/projects/mjolnir1/people/sjr729//tutorial/skimming_scripts/oedulis/fastq/skims_processing_pipeline/kraken/

GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/ncbi_dataset/data/GCF_947568905.1/GCF_947568905.1_xbOstEdul1.1_genomic.fna"

cd /projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/skmer2_with_skmer1lib_ref

python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/skmer_new-err.py --debug reference $INPUT_DIR -r $GENOME -p 2



# subsamp
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/ncbi_dataset/data/GCF_947568905.1/GCF_947568905.1_xbOstEdul1.1_genomic.fna"

# subsamp -10ind
RESPECT_SPECTRA=/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/respect-reference.txt
python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/skmer_new-err.py --debug reference /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/skims_processing_pipeline/subsamp -r $RESPECT_SPECTRA -p 2 -l ./skmer2_library_respect_2




# new prepocess no merging
module load parallel
conda activate skimming_echarvel
INPUT_DIR_2="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/skims_processing_pipeline_oct24/kraken/"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/ncbi_dataset/data/GCF_947568905.1/GCF_947568905.1_xbOstEdul1.1_genomic.fna"
cd /projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis
python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/skmer_new-err.py --debug reference $INPUT_DIR_2 -r $GENOME -l ./skmer2_library_ref_newpreprocess/ -p 2
