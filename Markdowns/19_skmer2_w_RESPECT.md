# Herring
# 1.SKMER1
```bash
INPUT_DIR=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/subsampled_4
OUTPUT_DIRECTORY=/projects/mjolnir1/people/sjr729/Skmer_ms/Herring/
conda activate skimming_echarvel
module load parallel

cd $OUTPUT_DIRECTORY
/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/run_skmer.sh -i $INPUT_DIR -t 20 -p 10 -o $OUTPUT_DIRECTORY/skmer1/skmer1_library
```

# 2.RESPECT
```bash
#it makes a directory with all the .hist files
#hsit_info file with the readlengths
#no need for run jellyfish

conda activate skimming_echarvel

NUM_THREADS=20
OUTPUT_DIRECTORY=/projects/mjolnir1/people/sjr729/Skmer_ms/Herring/
INPUT_DIR=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/subsampled_4

cd $OUTPUT_DIRECTORY

# Create separate directories for .hist files and hist_info.txt
mkdir --parents "${OUTPUT_DIRECTORY}/respect_4x/hist_files/"
echo -e "Input\tread_length" > "${OUTPUT_DIRECTORY}/respect_4x/hist_info.txt"

# Adjusted the path to match your directory structure
for directory in $(find "${OUTPUT_DIRECTORY}/skmer1/skmer1_library/skmer_library/" -maxdepth 1 -mindepth 1 -type d); do
    file=${directory##*/}
    
    # Ensure get_read_length function is available or replace with appropriate command
    read_len=$(grep 'read_length' "${directory}/${file}.dat" | awk '{print $2}')
    
    echo -e "${file}.hist\t${read_len}" >> "${OUTPUT_DIRECTORY}/respect_4x/hist_info.txt"
    ln --symbolic "$(realpath "${directory}/${file}.hist")" "${OUTPUT_DIRECTORY}/respect_4x/hist_files/"
done

# Now run respect, pointing to the hist_files directory
respect --input-directories "${OUTPUT_DIRECTORY}/respect_4x/hist_files/" \
        --info-file "${OUTPUT_DIRECTORY}/respect_4x/hist_info.txt" \
        --output-directory "${OUTPUT_DIRECTORY}/" \
        --spectra-output-size 50 \
        --iterations 1000 \
        --threads "${NUM_THREADS}"

```

# 3.SKMER2
```bash
cd /projects/mjolnir1/people/sjr729/Skmer_ms/Herring
#paste the estimated genome length to estimated-spectra.txt as a last column
paste <(cat estimated-spectra_1.txt) <(cut -f5 estimated-parameters_1.txt) > respect-reference.txt

RESPECT_SPECTRA=/projects/mjolnir1/people/sjr729/Skmer_ms/Herring/respect-reference.txt
INPUT_DIR=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/subsampled_4


python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/skmer_new-err.py --debug reference $INPUT_DIR -r $RESPECT_SPECTRA -l ./skmer2_library/ -p 2
```


# Oedulis
```bash
INPUT_DIR=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/skims_processing_pipeline/kraken/
OUTPUT_DIRECTORY=/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/
conda activate skimming_echarvel
module load parallel

cd $OUTPUT_DIRECTORY
/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/run_skmer.sh -i $INPUT_DIR -t 20 -p 20 -o $OUTPUT_DIRECTORY/skmer1/skmer1_library

```
# SKMER1
# RESPECT
# SKMER2
