## SRA to FastQ Conversion for HanClupea

Transforms .sra files in ERR* folders to fastq.gz files.

```bash
#!/bin/bash
base_dir="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testHanClupea"
for dir in "$base_dir"/SRR*/; do
    cd "$dir" || continue
    if ! ls *.fastq.gz 1> /dev/null 2>&1; then
        for file in ./*.sra; do
            fastq-dump --split-files "$file" && pigz "${file%.sra}"_*.fastq
        done
        rm -f ./*.sra
    else
        echo "Skipped: .fastq.gz files already exist in $dir"
    fi
    cd "$base_dir"
done
echo "All processing complete."
```

## Skmer Preprocessing Pipeline Launch for ClupeaAtmore

Executes Skmer preprocessing on fastq.gz files.

```bash
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/
conda activate tutorial
sbatch --wrap="bash ../skims_processing_pipeline.sh -x ./ -r 38 -f 38 > AtCluSkmin_sbatch_22jan.log"
```

```bash
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/
conda activate tutorial
bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea -r 39 -f 39 > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea_12feb24_screen.log"

```



## Fast Skmer Pipeline for Eduardo

Launches a faster Skmer preprocessing pipeline.

```bash
conda activate tutorial
bash ../fast_skims_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/ > Magpie_fasterskim_30jan.log
```

## File Movement Based on Matching Base Names for ClupeaAtmore

Moves .gz files to a target directory if their base names match.

```bash
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
for file in *.gz; do
    base_name=$(basename "$file" | cut -d '_' -f 1)
    if grep -q "${base_name}" ./skims_processing_pipeline/consult/*; then
        mv "$file" ./done
        echo "Moved $file to ./done"
    fi
done
```

## File Relocation from `./done` for Magpie

Relocates files if their first five characters do not match any in the target directory.

```bash
for f in ./done/*; do
    if ! ls ./skims_processing_pipeline/bbmap/${f:6:5}* 1> /dev/null 2>&1; then
        mv "$f" ./
    fi
done
```

## Skim Stats Refresh and Update

Concatenates .dat files to a .csv for postprocessing.

```bash
DATE=$(date +%d.%m)
cat ./library/*/*dat > "stats-postprocess_$DATE.csv"
```

## Preprocess Stats and Output Wrangling for ClupeaAtmore

Generates stats and identifies bad specimens.

```bash
DATE=$(date +%d.%m)
bash ../post_processing_pipeline.sh ../ testClupea postprocess_18jan library
python csv_to_tsv_stats.py "stats-postprocess_$DATE.csv" "stats-postprocess_$DATE.tsv"
python badspecimens_identifier.py "stats-postprocess_$DATE.tsv" library library_badindividuals
```

## Skmer Distance Calculation

Computes distance matrices using Skmer.

```bash
DATE=$(date +%d.%m)
skmer distance library -t -o "jc-$DATE-dist-mat"
```

## Subsample to 4x and estimate with Skmer
```bash
cd /path/to/dir
conda activate tutorial 
module load parallel
bash subsample_and_estimate.sh -i skims_processing_pipeline/kraken -c 4 -t 40
```

## Quick Phylogeny with FastTree

Performs phylogenetic analysis using FastTree.

```bash
module load fasttree/2.1.11
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
bash ../tsv_to_phymat.sh "jc-$DATE-dist-mat.txt" "jc-$DATE-dist-mat.phy"
../fastme-2.1.5/binaries/fastme-2.1.5-linux64 -i "jc-$DATE-dist-mat.phy" -o "backbone-fastme_jc-$DATE-dist-mat.tre"
../fastme-2.1.5/binaries/fastme-2.1.5-linux64 -i "jc-$DATE-dist-mat.phy" -u "backbone-fastme_jc-$DATE-dist-mat.tre" -o "backbone-fastme_jc-$DATE-dist-mat.tre"
nw_reroot "backbone-fastme_jc-$DATE-dist-mat.tre" | nw_display -
```
