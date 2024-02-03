## The scripts transform .sra files in ERR* folder into fastq.gz files in 
# For HanClupea

```bash
!/bin/bash

# Directory containing the folders with .sra files
base_dir="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testHanClupea"

# Loop through each directory starting with SRR*
for dir in "$base_dir"/SRR*/; do
    echo "Checking directory: $dir"

    # Move into the directory, or skip if it doesn't exist
    cd "$dir" || continue

    # Check for existing .fastq.gz files in the directory
    if ls *.fastq.gz 1> /dev/null 2>&1; then
        echo "Skipped: .fastq.gz files already exist in $dir"
        # Return to the base directory and continue with next directory
        cd "$base_dir" || exit
        continue
    fi

    echo "Processing directory: $dir"

    # Enable globbing to match file names
    shopt -s nullglob
    shopt -s dotglob

    # Find .sra files in the directory
    sra_files=(./*.sra)

    if [ ${#sra_files[@]} -gt 0 ]; then
        # Convert .sra files to FASTQ and compress to .gz format
        for file in "${sra_files[@]}"; do
            accession=$(basename "$file" .sra)
            fastq-dump --split-files "$file" && pigz "${accession}"_*.fastq
        done

        # Remove .sra files after conversion
        rm -f ./*.sra
    fi

    # Return to the base directory
    cd "$base_dir" || exit
done

echo "All processing complete."
```


## The script is launching Skmer preprocessing pipeline on /path/to/folder/*.fastq.gz
# For ClupeaAtmore

```bash
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/
conda activate tutorial
sbatch --job-name=AtCluSkmin_sbatch --output=AtCluSkmin_sbatch_22jan.out --error=AtCluSkmin_sbatch_22jan.err --ntasks=1 --cpus-per-task=40 --mem=180G --time=80:00:00 --mail-type=begin --mail-type=end --mail-type=fail --mail-user=homerejalves.monteiro@sund.ku.dk --wrap="cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea && bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/ -r 38 -f 38 >/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/AtCluSkmin_sbatch_22jan.log "


sbatch --job-name=HanCluSkmin_sbatch --output=HanCluSkmin_sbatch.out --error=HanCluSkmin_sbatch.err --ntasks=1 --cpus-per-task=40 --mem=180G --time=80:00:00 --mail-type=begin --mail-type=end --mail-type=fail --mail-user=homerejalves.monteiro@sund.ku.dk --wrap="cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testHanClupea && bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testHanClupea/ -r 38 -f 38 >/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testHanClupea/AtCluSkmin_sbatch_22jan.log "

```

## The script is launching Skmer faster **Eduardo** preprocessing pipeline on /path/to/folder/*.fastq.gz

 ```bash
 conda activate tutorial

 bash ../fast_skims_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/ -r 38 -f 38 >/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/Magpie_fasterskim_30jan.log
```


## The script moves .gz files from a source directory to a target directory if their base names (derived from the filenames before the first underscore) match with any files present in the target directory.
# For ClupeaAtmore

```bash
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea

# Directory containing the files to be moved
SOURCE_DIR="./"

# Directory containing the files to match against
TARGET_DIR="./skims_processing_pipeline/consult/"

# Loop over .gz files in the source directory
for file in *.gz; do base_name=$(basename "$file" | cut -d '_' -f 1); if ls -l ./skims_processing_pipeline/consult/ | grep -q "${base_name}"; then mv "$file" ./done; echo "Moved $file to ./done"; fi; done
```

## Moves files from the `./done` directory to the current directory (`./`) if their first five characters do not match the first five characters of any file name in the `./skims_processing_pipeline/bbmap/` directory.
# for magpie
```bash
for f in ./done/*; do if ! ls ./skims_processing_pipeline/bbmap/${f:6:5}* 1> /dev/null 2>&1; then mv "$f" ./; fi; done
```

## Refresh and update the skim stats .csv 
```bash
DATE=$(date +%d.%m)

cat ./library/*/*dat > stats-postprocess_$DATE.csv
```


## Get stats of the preprocess and wrangling the output with the two python scripts
# For ClupeaAtmore
```bash
DATE=$(date +%d.%m)

bash ../post_processing_pipeline.sh ../ /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea postprocess_18jan /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/library

module load python

python csv_to_tsv_stats.py stats-postprocess_$DATE.csv stats-postprocess_$DATE.tsv

python badspecimens_identifier.py stats-postprocess_$DATE.tsv library library_badindividuals
```




## The script is launching Skmer distance to get dist-mat and jc-dist-mat (Jukes-Cantor correction)
```bash
DATE=$(date +%d.%m)
skmer distance library -t -o jc-$DATE-dist-mat
```

## The script is doing a quick phylogeny using fasttree
```bash
module load fasttree/2.1.11
#On a subset for Clupea dataset
#Skmer distance on existing library
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
bash ../tsv_to_phymat.sh jc-$DATE-dist-mat jc-$DATE-dist-mat.phy
# Infer backbone from scratch
../fastme-2.1.5/binaries/fastme-2.1.5-linux64  -i jc-$DATE-dist-mat.phy  -o backbone-fastme_jc-$DATE-dist-mat.tre
# Recompute branch lengths:
../fastme-2.1.5/binaries/fastme-2.1.5-linux64  -i jc-$DATE-dist-mat.phy  -o backbone-fastme_jc-$DATE-dist-mat.tre -u backbone-fastme_jc-$DATE-dist-mat.tre
#
nw_reroot backbone-fastme_jc-$DATE-dist-mat.tre |nw_display -
```
