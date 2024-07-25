## SRA to FastQ Conversion 
# For oedulis
# Transforms .sra files in ERR* folders to fastq.gz files.
```bash
#!/bin/bash
base_dir="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis"
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

# Alternative version with.csv file containing run accession in the first column Stein Mortensen collected oedulis
```bash
module load sratoolkit sra-tools pigz 
base_dir="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq"
# Path to the CSV file
csv_file="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/oedulis_sra.csv"
# Navigate to the base directory
cd "$base_dir" || exit
# Read SRA run IDs sampled by Stein Mortensen from the CSV, skip the header line
sra_ids=$(tail -n +2 "$csv_file" | cut -d',' -f1,11 | grep "Stein Mortensen" | cut -d',' -f1)
#sra_ids=$(tail -n +2 "$csv_file" | cut -d',' -f1,11 | grep "Ane" | cut -d',' -f1)

for sra_id in $sra_ids; do
    # Check if .fastq.gz files already exist for this SRA ID
    if [ -f "${sra_id}_1.fastq.gz" ] && [ -f "${sra_id}_2.fastq.gz" ]; then
        echo "Skipped: .fastq.gz files already exist for $sra_id"
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
        else
            echo "Error: SRA file for $sra_id not found."
        fi
    fi
done
echo "All processing complete."
```




## Stein Mortensen sampled skim_processed .fq in fastq/
```bash
base_dir="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/Stein_Mortensen_fqs/"
# Path to the CSV file
csv_file="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/oedulis_sra.csv"

# Navigate to the base directory
cd "$base_dir" || exit

# Read SRA run IDs sampled by Stein Mortensen from the CSV, skip the header line
sra_ids=$(tail -n +2 "$csv_file" | grep "Stein Mortensen" | cut -d',' -f1)

for sra_id in $sra_ids; do
    # Check if .fq files already exist for this SRA ID in oedulis/fastq/Stein_Mortensen_fqs/ if not, cp the file in oedulis/fastq/Stein_Mortensen_fqs/
    if [ -f "unclassified-kra_${sra_id}_.fq" ] ; then
        echo "Skipped: .fq files already exist for $sra_id"
    else
        # Check if the SRA file was downloaded and exists
        if [ -f "../skims_processing_pipeline/kraken/unclassified-kra_${sra_id}_.fq" ]; then
            cp "../skims_processing_pipeline/kraken/unclassified-kra_${sra_id}_.fq" ./
        else
            echo "File not found: unclassified-kra_${sra_id}_.fq"
        fi
    fi
done
```

## Skmer Preprocessing Pipeline 
# Launch for Herring
Executes Skmer preprocessing on fastq.gz files.

```bash
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/
conda activate Mar_skmer_pip
sbatch --wrap="bash ../skims_processing_pipeline.sh -x ./ -r 38 -f 38 > AtCluSkmin_sbatch_22jan.log"

cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/
conda activate Mar_skmer_pip
module load parallel kraken2 respect consult-ii
bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/done -r 39 -f 39 > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/log/27may24_screen_skimprocess.log


cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/
conda activate Mar_skmer_pip
module load parallel kraken2 respect consult-ii
bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/rawfastqmodernClupea/ -r 40 -f 40 > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea_9apr24_screen2.log
```

# For oedulis
```bash
conda activate skmer_2_test
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/
module load parallel kraken2 respect consult-ii
bash ../../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq -r 40 -f 40 > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/oedu_preprocess_29jun24_screen2.log
```


## Fast Skmer Pipeline for Eduardo

# Launches a faster Skmer preprocessing pipeline.

```bash
conda activate tutorial
bash ../fast_skims_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/ > Magpie_fasterskim_9apr24.log
```


## Wrangling on skim_preprocessing output

# File Movement Based on Matching Base Names for ClupeaAtmore
# Moves .gz files to a target directory if their base names match.
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

# File Relocation from `./done` for Magpie
# Relocates files if their first five characters do not match any in the target directory.
```bash
for f in ./done/*; do
    if ! ls ./skims_processing_pipeline/bbmap/${f:6:5}* 1> /dev/null 2>&1; then
        mv "$f" ./
    fi
done
```

# Concatenates .dat files to a .csv for postprocessing.
```bash
DATE=$(date +%d.%m)
cat ./library/*/*dat > "stats-postprocess_$DATE.csv"
```
# specific for Clupea skmer raw cov
```bash
for file in ./library/done/*/*.dat; do
    awk -v fname=$(basename "$file" .dat) '{print fname"    "$0}' "$file"
done > stats-rowcov_postprocess_$DATE.csv
```
# specific for Magpie skmer raw cov
```bash
for file in ./library/*/*.dat; do
    awk -v fname=$(basename "$file" .dat) '{print fname"    "$0}' "$file"
done > stats-rowcov_postprocess_$DATE.csv
```


## Preprocess Stats and Output Wrangling for ClupeaAtmore

# Generates stats and identifies bad specimens.

```bash
DATE=$(date +%d.%m)
bash ../post_processing_pipeline.sh ../ testClupea postprocess_18jan library
python csv_to_tsv_stats.py "stats-postprocess_$DATE.csv" "stats-postprocess_$DATE.tsv"
python badspecimens_identifier.py "stats-postprocess_$DATE.tsv" library library_badindividuals
```





## Skmer Distance Calculation with skmer1

# Computes distance matrices using Skmer.
```bash
conda activate skmer_2_test
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie
module load parallel
DATE=$(date +%d.%m)
skmer reference /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/skims_processing_pipeline/kraken/ -t -o "jc-$DATE-dist-mat" -p 29 -o "skmer1_Ref-$DATE-dist-mat" -l library_skmer1_Ref_$DATE
```


## Skmer2 Distance Calculation with skmer2 for Magpie

```bash
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie
conda activate skmer_2_test
conda install jellyfish seqtk mash
module load parallel
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/ncbi_dataset/data/GCA_013398635.1/GCA_013398635.1_ASM1339863v1_genomic.fna"
DATE=$(date +%d.%m)
# build the skmer2 library
python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/__main_TESTING.py --debug reference /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/skims_processing_pipeline/kraken/ -r $GENOME -p 2 -o "skmer2_Magpie_Ref-$DATE-dist-mat" -l library_skmer2_Magpie_Ref

```


## subsample_and_estimate.sh to 4x and estimate with Skmer
```bash
#cd /path/to/dir
conda activate Mar_skmer_pip 
module load parallel kraken2 
DATE=$(date +%d.%m)
#bash subsample_and_estimate.sh -i skims_processing_pipeline/kraken -c 4 -t 40
bash subsample_and_estimate.sh -i skims_processing_pipeline/kraken -c 2 -t 40 #not complete for Magpie
bash subsample_and_estimate.sh -i skims_processing_pipeline/kraken -c 1 -t 40 #done
bash subsample_and_estimate.sh -i skims_processing_pipeline/kraken -c 0.5 -t 40  > .log #done
bash subsample_and_estimate.sh -i skims_processing_pipeline/kraken -c 0.5 -t 40  > .log #done
```

## Troubleshoot subsample_and_estimate.sh
```bash
# For testMagpie
bash subsample_and_estimate.sh -i skims_processing_pipeline/kraken -c 1 -t 40 -o OUT_subsample_and_estimate_1x  > "$DATE"_test_toubleshoot_11apr_1x.log
bash subsample_and_estimate.sh -i skims_processing_pipeline/kraken -c 1 -t 40 -o OUT_subsample_and_estimate_1x  > "$DATE"_test_toubleshoot_withmodules_1x.log
bash subsample_and_estimate.sh -i skims_processing_pipeline/kraken -c 0.5 -t 40 -o OUT_subsample_and_estimate_0.5x  > "$DATE"_test_toubleshoot_withmodules_0.5x.log

#For testClupea
conda activate Mar_skmer_pip 
module load parallel kraken2 
DATE=$(date +%d.%m)
#bash ../testMagpie/subsample_and_estimate.sh -i skims_processing_pipeline_jan24/kraken -c 1 -t 40 -o OUT_subsample_and_estimate_1x  > "$DATE"_testClupea_toubleshoot_withmodules_1x.log #### subsample 1x done
bash subsample_and_estimate.sh -i skims_processing_pipeline_jan24/kraken -c 0.25 -t 40 -o OUT_subsample_and_estimate_0.25x  > "$DATE"_testClupea_toubleshoot_withmodules_0.25x.log

for file in skims_processing_pipeline_jan24/kraken/*__merged; do
  mv "$file" "${file}.fq"
done

# I have cp the updated subsample_and_estimate.sh in testClupea 24apr24
bash ./subsample_and_estimate.sh -i ./skims_processing_pipeline_jan24/kraken -t 30 -o ./OUT_subsample_and_estimate_1x -c 1 2>&1 > "$DATE"_testClupea_subsample_n_estimate_troubleshoot.log


sbatch --job-name=ClupeaEstimate_4x \
--output=Estimate_4x_%j.out \
--error=Estimate_4x_%j.err \
--ntasks=1 \
--cpus-per-task=40 \
--mem=180G \
--time=45:00:00 \
--mail-type=BEGIN,END,FAIL \
--mail-user=homerejalves.monteiro@sund.ku.dk \
--wrap="conda activate Mar_skmer_pip; module load parallel kraken2; DATE=\$(date +%d.%m); bash /subsample_and_estimate.sh -i skims_processing_pipeline_jan24/kraken -c 2 -t 40 -o OUT_subsample_and_estimate_2x > \${DATE}_testClupea_subsample_sbatch_2x.log"



#template
bash ./template_subsample_27apr.sh -i ./skims_processing_pipeline_jan24/kraken -t 30 -o ./OUT_subsample_and_estimate_4x_template -c 4 2>&1 > "$DATE"_testClupea_templatesubsample_.log



#30apr24
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
bash ./subsample_and_estimate.sh -i ./skims_processing_pipeline_jan24/kraken -t 30 -o ./OUT_subsample_and_estimate_4x_template -c 4 2>&1 > "$DATE"_testClupea_templatesubsample.log


#1may
conda activate skmer_2_test
DATE=$(date +%d.%m)
module load parallel kraken2 
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
bash ./subsample_and_estimate.sh -i ./skims_processing_pipeline_jan24/kraken -t 30 -o ./OUT_subsample_and_estimate_4x_template -c 4 2>&1 > "$DATE"_skmer2_Clupea_templatesubsample.log
bash ./subsample_and_estimate.sh -i ./skims_processing_pipeline_jan24/kraken -t 18 -o ./OUT_subsample_and_estimate_2x_skmer2 -c 2 2>&1 > "$DATE"_skmer2__Clupea_2xubsample.log


#2may
conda activate Mar_skmer_pip 
DATE=$(date +%d.%m)
module load parallel kraken2 
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
bash ./subsample_and_estimate.sh -i ./skims_processing_pipeline_jan24/kraken -t 18 -o ./OUT_subsample_and_estimate__Mar_skmer_pip_2x_skmer2 -c 2 2>&1 > "$DATE"_skmer2_condaMar_skmer_pip-_Clupea_2xubsample.log


#3may
conda activate Mar_skmer_pip 
DATE=$(date +%d.%m)
module load parallel kraken2 
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
bash ./subsample_and_estimate.sh -i ./skims_processing_pipeline_jan24/kraken -t 18 -o ./OUT_subsample_and_estimate__Mar_skmer_pip_0.25x_skmer2 -c 0.25 2>&1 > "$DATE"_skmer2_condaMar_skmer_pip-_Clupea_0.25xubsample.log


du --summarize --block-size=1K ./skims_processing_pipeline_jan24/kraken/* | awk '$1 <= 10485760' | cut -f 2 | xargs

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


## sbatch-ready script for preprocessing pipeline with skmer
```bash
#With conda env= tutorial 
sbatch --job-name=4xSkmin_sbatch --output=4xSkmin_sbatch_16jan.out --error=4xSkmin_sbatch_16jan.err --ntasks=1 --cpus-per-task=40 --mem=180G --time=400:00:00 --mail-type=begin --mail-type=end --mail-type=fail --mail-user=homerejalves.monteiro@sund.ku.dk --wrap="cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea && conda activate tutorial && bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea -r 40 -f 40 "

#With conda env= tutorial 
With conda env= Mar_skmer_pip 
sbatch --job-name=Skmin_sbatch_27mar --output=Skmin_sbatch_27mar.out --error=Skmin_sbatch_27mar.err --ntasks=1 --cpus-per-task=40 --mem=180G --time=300:00:00 --mail-type=begin --mail-type=end --mail-type=fail --mail-user=homerejalves.monteiro@sund.ku.dk --wrap="cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea && conda activate Mar_skmer_pip && bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea -r 40 -f 40 "

#Alt
sbatch --job-name=Skmin_sbatch_27mar --output=Skmin_sbatch_27mar.out --error=Skmin_sbatch_27mar.err --ntasks=1 --cpus-per-task=40 --mem=180G --time=300:00:00 --mail-type=begin --mail-type=end --mail-type=fail --mail-user=homerejalves.monteiro@sund.ku.dk --wrap="cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea && source activate Mar_skmer_pip && bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea -r 40 -f 40"

```


## Kraken only script after consult decontamination for Clupea dataset
```bash
#!/bin/bash

# Define the directories
input_dir="skims_processing_pipeline_jan24/consult/"  # Change this to your actual consult directory path
kraken_db="../../skimming_scripts/kraken2/krakenlib/"          # Change this to your Kraken2 database path
output_dir="skims_processing_pipeline_jan24/kraken"     # Specify your desired output directory for Kraken results

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop over each file in the consult directory
for file in "$input_dir"/*; do
  filename=$(basename -- "$file")
  output_file="${output_dir}/unclassified-kra_${filename}"

  # Run Kraken2
  echo "Running Kraken2 decontamination on $filename"
  kraken2 --db "$kraken_db" "$file" --unclassified-out "$output_file"
  echo "Kraken2 decontamination completed for $filename"
done

echo "All Kraken2 decontamination processes are complete."

```

## Skmer operations, Respect pipeline, and Skmer distance calculation after consult + kraken decontamination
```bash
#!/bin/bash
module load parallel kraken2 skmer respect

# Paths
input_dir="/path/to/skims_processing_pipeline/kraken"
lib_dir="/path/to/library"
output_dir="/path/to/output"
temp_dir="/path/to/temp"

# Settings
cores=8
iterations=1000
threads=8

# Ensure necessary directories exist
mkdir -p "$lib_dir" "$output_dir" "$temp_dir"

# Skmer Operations + Respect Pipeline
for file in "$input_dir"/*_read.fq; do  #*_read.fq -> this might not be the suffix
    filename=$(basename "$file")
    sample_id="${filename%_read.fq}"

    # Skmer Operations
    echo "Running Skmer operations for $sample_id"
    if [ ! -d "${lib_dir}/${sample_id}" ]; then
        mkdir -p "${lib_dir}/${sample_id}"
        # Placeholder for Skmer sketching command, replace with actual Skmer command
        skmer sketch "$file" -o "${lib_dir}/${sample_id}" --threads $cores
    fi

    # Respect Pipeline (Pseudo-code, replace with actual commands if applicable)
    echo "Running Respect pipeline for $sample_id"
    if [ ! -d "${output_dir}/respect/${sample_id}" ]; then
        mkdir -p "${output_dir}/respect/${sample_id}"
        # Placeholder for Respect command, replace with actual Respect command
        respect -i "$file" -l "${lib_dir}/${sample_id}" -o "${output_dir}/respect/${sample_id}" --iterations $iterations --threads $threads
    fi

    echo "Completed processing for $sample_id"
done

# Skmer Distance Calculation
echo "Calculating Skmer distances"
# Assuming skmer distance calculation requires a list of all sketches directories
skmer_lib_dirs=$(find "$lib_dir" -type d -name 'unclassified-kra_*' -print)
skmer distance $skmer_lib_dirs --threads $cores -o "$output_dir/skmer_distances"

# Post-Processing
echo "Starting post-processing steps"
# Assuming a generic post-processing command, replace with actual needs
../post_processing_pipeline.sh "$output_dir/skmer_distances" -o "$output_dir/final_results"

echo "Cleaning up temporary files"
rm -rf "$temp_dir"

echo "All operations completed successfully"
```


## Skmer1 with script run_skmer.sh to get distance matrix of skim_preprocessed reads
# oedulis
```bash
conda activate skmer_2_test
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis
module load parallel
DATE=$(date +%d.%m)
/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/run_skmer.sh -i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/skims_processing_pipeline/kraken/ -t 19 -p 10 -o /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/skmer/library_runskmer1_oedulis_condaskmer2test_$DATE
 ```
# Magpie
```bash 
conda activate skmer_2_test
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie
module load parallel
DATE=$(date +%d.%m)
/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/run_skmer.sh -i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/skims_processing_pipeline/kraken/ -t 40 -p 20 -o library_runskmer1_magpie_condaskmer2test_$DATE
```
# Herring
```bash
conda activate skmer_2_test
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
module load parallel
DATE=$(date +%d.%m)
/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/run_skmer.sh -i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skims_processing_pipeline_jan24/kraken -t 29 -p 10 -o library_clupea_run_skmer1_Ref_$DATE
```

## Use run_skmer.sh to get distance matrix of skimmed preprocessed read mapped to the ref genome
# Herring
```bash
conda activate skmer_2_test
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
module load parallel
DATE=$(date +%d.%m)
 /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/run_skmer.sh -i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/mapped_fq -t 19 -p 10 -o library_runskmer_mappedreads_skmer1_condaskmer2test_$DATE
 # result are in library_runskmer_mappedreads_skmer1_condaskmer2test_26.06

```

# For subsampled at 4x herring 
```bash
conda activate skmer_2_test
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
module load parallel
DATE=$(date +%d.%m)
/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/run_skmer.sh -i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/subsampled_4/ -t 29 -p 10 -o library_subsampled_4_skmer1_Ref_$DATE
```


## Skmer2 new scripts Jun24 on herring skim_preprocessed + mapped reads  with genome ref
```bash
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
conda activate skmer_2_test
module load parallel
DATE=$(date +%d.%m)
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/skmer_new-err.py --debug reference /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/mapped_fq/ -r $GENOME -p 4 -o "skmer2_echarvel_Clupea_Ref-$DATE-dist-mat" -l library_skmer2_echarvel_Clupea_Ref-$DATE
```

## subsample_and_estimate2  (needs module parallel to work)
# Magpie skim_preprocessed reads
```bash
conda activate skmer_2_test 
DATE=$(date +%d.%m)
module load parallel
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie
INPUTDIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/skims_processing_pipeline/kraken/"
bash /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/subsample_and_estimate2.sh -i $INPUTDIR -t 16 c 4,2,1,0.5,0.25 2>&1 > subsample2"$DATE"_condaskmer2test-_Magpieskim_preprocessed.log


```


## skmer2 with skmer_new-err.py with -p 1 on Herring skims preprocessed + mapped reads to ref Genome 
```bash
conda activate skmer_2_test
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
module load parallel
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
DATE=$(date +%d.%m)
python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/skmer_new-err.py --debug reference /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/mapped_fq/ -r $GENOME -p 1 -o "skmer2_echarvel_Clupea_Ref-$DATE-dist-mat" -l library_skmer2_echarvel_Clupea_Ref-$DATE
```