# Mapping reads to GENOME SANDEEL
```bash

# Create dependencies for the genome
samtools faidx $GENOME
bwa index $GENOME
picard CreateSequenceDictionary R=$GENOME O=${GENOME%.fna}.dict


# Load necessary modules
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18 jdk/1.8.0_291 picard/2.27.5 parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15 python/3.11.3 openjdk/17.0.8 gatk/3.8 bedtools bbmap mosdepth

# Define directories and genome reference
DATAINPUT="/projects/diversea/data/20240809_Glenn_sand_eels_data/HN00220216_hdd1/RawData/"
DATAOUTPUT="/projects/diversea/people/sjr729/Sandeels/mapped/RawReadsMappedToGenome"
GENOME="/projects/diversea/people/sjr729/Sandeels/Genome/GCA_026929865.2_TBG_H_lanceolatus_asm_v1.1_genomic.fna"
BBMAP_DIR="/projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/bbmap"
FASTQ_OUTPUT="/projects/diversea/people/sjr729/Sandeels/mapped/RawReadsMappedToGenome/fastq"
BBMAP_OUTPUT="/projects/diversea/people/sjr729/Sandeels/mapped/RawReadsMappedToGenome/bbmap_fastq"
STATS_OUTPUT="/projects/diversea/people/sjr729/Sandeels/mapped/RawReadsMappedToGenome/individual_stats"
GLOBAL_STATS_OUTPUT="/projects/diversea/people/sjr729/Sandeels/mapped/RawReadsMappedToGenome/global_mapping_stats.tsv"

# Create necessary output directories
mkdir -p "$DATAOUTPUT" "$FASTQ_OUTPUT" "$BBMAP_OUTPUT" "$STATS_OUTPUT"
cd $DATAOUTPUT

# Initialize the global stats file with a header if it doesn't exist
if [ ! -f "$GLOBAL_STATS_OUTPUT" ]; then
    echo -e "Sample\tRaw_Reads\tRaw_Bases\tTrimmed_Bases\tMapped_Reads_BWA\tPercentage_Mapped_Reads_BWA\tMapped_Reads_Post_BWA_BBMAP\tMapped_Percentage_Post_BWA_BBMAP\tMosdepth_Avg_Depth\tProp_g_cov_1\tProp_g_cov_2\tProp_g_cov_3\tProp_g_cov_4" > "$GLOBAL_STATS_OUTPUT"
fi


# Memory and split parameters for bbmap
dukmem=16
splitsize=20000000  # Number of lines per chunk

# Loop through each _1.fastq.gz file in the DATAINPUT directory
for filepath1 in "$DATAINPUT"/*_1.fastq.gz; do
    filepath2="${filepath1/_1.fastq.gz/_2.fastq.gz}"
    filename=$(basename "$filepath1")
    base=${filename%_1.fastq.gz}

    # Overwrite individual stats file at the start of each sample iteration
    individual_stats_file="$STATS_OUTPUT/${base}_stats.tsv"
    echo -e "Sample\tRaw_Reads\tRaw_Bases\tTrimmed_Bases\tMapped_Reads_BWA\tPercentage_Mapped_Reads_BWA\tMapped_Reads_Post_BWA_BBMAP\tMapped_Percentage_Post_BWA_BBMAP\tMosdepth_Avg_Depth\tProp_g_cov_1\tProp_g_cov_2\tProp_g_cov_3\tProp_g_cov_4" > "$individual_stats_file"


    ##############################
    # Step 1: Calculate raw reads and raw bases before mapping
    ##############################
    echo "Calculating raw reads and raw bases for $base"

    # Calculate raw reads
    raw_reads_forward=$(zcat "$filepath1" | wc -l | awk '{print $1/4}')
    raw_reads_reverse=$(zcat "$filepath2" | wc -l | awk '{print $1/4}')
    raw_reads=$((raw_reads_forward + raw_reads_reverse))
    
    # Calculate raw bases
    raw_bases_forward=$(zcat "$filepath1" | awk 'NR%4==2' | tr -d "\n" | wc -m)
    raw_bases_reverse=$(zcat "$filepath2" | awk 'NR%4==2' | tr -d "\n" | wc -m)
    raw_bases=$((raw_bases_forward + raw_bases_reverse))

    # Initialize trimmed bases as 0 (for now, can update if trimming step is added later)
    trimmed_bases=0

    # Add raw reads and raw bases to the individual stats file
    echo -e "$base\t$raw_reads\t$raw_bases\t$trimmed_bases" >> "$individual_stats_file"

    ##############################
    # Step 2: Mapping raw reads to the GENOME using bwa mem
    ##############################
    if [ -f "$DATAOUTPUT/$base.sort.minq20.bam" ]; then
        echo "Skipping mapping for $base, BAM file already exists."
    else
        echo "Mapping raw reads to GENOME for $base"
        bwa mem -t 30 -R "@RG\tID:$base\tSM:$base\tPL:Illumina" "$GENOME" "$filepath1" "$filepath2" > "$DATAOUTPUT/$base.sam" 2> "$DATAOUTPUT/$base.bwa.log"
        
        if [ $? -ne 0 ]; then
            echo "BWA failed for $base"
            continue
        fi

        samtools view -bS -h -q 20 -F 4 "$DATAOUTPUT/$base.sam" > "$DATAOUTPUT/$base.bam"
        samtools sort "$DATAOUTPUT/$base.bam" -o "$DATAOUTPUT/$base.sort.minq20.bam"
        samtools index "$DATAOUTPUT/$base.sort.minq20.bam"
        rm "$DATAOUTPUT/$base.sam" "$DATAOUTPUT/$base.bam"

        if [ $? -ne 0 ]; then
            echo "SAM/BAM processing failed for $base"
            continue
        fi
    fi

    ##############################
    # Step 3: BAM to Fastq conversion using bedtools bamtofastq
    ##############################
    fq_file="$FASTQ_OUTPUT/$base.fq"
    if [ -f "$fq_file" ]; then
        echo "Skipping BAM to Fastq conversion for $base, Fastq file already exists."
    else
        echo "Converting BAM to Fastq for $base"
        bedtools bamtofastq -i "$DATAOUTPUT/$base.sort.minq20.bam" -fq "$fq_file"
        if [ $? -ne 0 ]; then
            echo "BAM to Fastq conversion failed for $base"
            continue
        fi
    fi

    ##############################
    # Step 4: Process Fastq with bbduk and dedupe
    ##############################
    dedupe_merged_file="${BBMAP_OUTPUT}/${base}_bbmap.fq"
    dedupe_bam_file="${DATAOUTPUT}/${base}_bbmap_sorted.bam"  # BAM after bbduk and dedupe

    # Check if the deduplicated BAM file exists, and if not, create it
    if [ -f "$dedupe_merged_file" ] && [ ! -f "$dedupe_bam_file" ]; then
        echo "bbmap Fastq file exists, but BAM file is missing. Creating BAM file."
        bwa mem -t 30 "$GENOME" "$dedupe_merged_file" | samtools view -bS - | samtools sort -o "$dedupe_bam_file"
        samtools index "$dedupe_bam_file"
    elif [ -f "$dedupe_bam_file" ]; then
        echo "Skipping bbmap processing for $base, bbmap processed file and BAM file already exist."
    else
        # If neither the Fastq nor the BAM file exists, process everything
        echo "Processing Fastq with bbmap for $base"
        bbduk_output_prefix="${BBMAP_OUTPUT}/${base}_chunk"
        split -l "$splitsize" -d "$fq_file" "${bbduk_output_prefix}_"

        # Processing chunks with bbduk and dedupe
        for chunk in "${bbduk_output_prefix}_"*; do
            bbduk_out="${chunk}_bbduk.fq"
            dedupe_out="${chunk}_dedupe.fq"

            # Run bbduk on each chunk
            $BBMAP_DIR/bbduk.sh -Xmx${dukmem}g -Xms${dukmem}g in="$chunk" out="$bbduk_out" ref=adapters,phix ktrim=r k=23 mink=11 hdist=1 tpe tbo overwrite=true
            
            if [ -s "$bbduk_out" ]; then
                echo "bbduk completed successfully for $chunk."
            else
                echo "bbduk failed for $chunk."
                continue
            fi

            # Run dedupe on the bbduk output
            $BBMAP_DIR/dedupe.sh -Xmx${dukmem}g -Xms${dukmem}g in="$bbduk_out" out="$dedupe_out" overwrite=true
            
            if [ -s "$dedupe_out" ]; then
                echo "dedupe completed successfully for $chunk."
            else
                echo "dedupe failed for $chunk."
                continue
            fi

            # Clean up intermediate files (optional)
            rm "$bbduk_out" "$chunk"
        done

        # Merge all deduplicated chunks into one output file
        cat "${bbduk_output_prefix}_"*_dedupe.fq > "$dedupe_merged_file"
        rm "${bbduk_output_prefix}_"*_dedupe.fq

        # Convert the merged deduplicated Fastq to BAM
        bwa mem -t 30 "$GENOME" "$dedupe_merged_file" | samtools view -bS - | samtools sort -o "$dedupe_bam_file"
        samtools index "$dedupe_bam_file"

        echo "bbmap processing completed for $base."
    fi

    ##############################
    # Step 5: Calculate mapped reads after bwa mem using samtools flagstat
    ##############################
    mapped_reads_bwa=$(samtools flagstat "$DATAOUTPUT/$base.sort.minq20.bam" | grep "mapped (" | awk '{print $1}')
    percentage_mapped_reads_bwa=$(awk "BEGIN {printf \"%.2f\", ($mapped_reads_bwa * 100) / $raw_reads}")

    # Append mapped reads and percentage post-bwa to individual stats file
    sed -i "/^$base/s/$/\t$mapped_reads_bwa\t$percentage_mapped_reads_bwa/" "$individual_stats_file"

    ##############################
    # Step 6: Calculate mapped reads after bbduk and dedupe (Post-BWA BBMAP Stats)
    ##############################
    if [ -f "$dedupe_bam_file" ]; then
        post_bwa_bbmap_mapped_reads=$(samtools flagstat "$dedupe_bam_file" | grep "mapped (" | awk '{print $1}')
        post_bwa_bbmap_mapping_percentage=$(awk "BEGIN {printf \"%.2f\", ($post_bwa_bbmap_mapped_reads * 100) / $raw_reads}")


        # Append the stats for post-bwa bbmap mapping
        sed -i "/^$base/s/$/\t$post_bwa_bbmap_mapped_reads\t$post_bwa_bbmap_mapping_percentage/" "$individual_stats_file"
    else
        echo "No BAM file after bbmap processing found for $base"
        continue
    fi

    ##############################
    # Step 7: Calculate depth of coverage using mosdepth
    ##############################
    depth_summary="${DATAOUTPUT}/${base}.mosdepth.mosdepth.summary.txt"
    if grep -q "$base" "$GLOBAL_STATS_OUTPUT"; then
        echo "Skipping mosdepth calculation for $base, stats already exist."
    else
        mosdepth --threads 2 --by 10000 "$DATAOUTPUT/$base.mosdepth" "$dedupe_bam_file"
        average_depth=$(awk 'NR==2 {print $4}' "$depth_summary")

        # Append average depth to the stats files
        sed -i "/^$base/s/$/\t$average_depth/" "$individual_stats_file"
    fi


    ##############################
    # Step 8: Calculate proportion of genome covered by at least 1, 2, 3, and 4 reads using grep and awk
    ##############################
    base_coverage_file="$DATAOUTPUT/${base}_coverage.txt"

    if [ ! -f "$base_coverage_file" ]; then
        bedtools genomecov -ibam "$dedupe_bam_file" -g "$GENOME.fai" > "$base_coverage_file"
    fi

    prop_gcov_1=$(cat "$base_coverage_file" | grep "genome" | awk '$2 >= 1 {sum += $5} END {printf "%.2f\n", sum * 100}')
    prop_gcov_2=$(cat "$base_coverage_file" | grep "genome" | awk '$2 >= 2 {sum += $5} END {printf "%.2f\n", sum * 100}')
    prop_gcov_3=$(cat "$base_coverage_file" | grep "genome" | awk '$2 >= 3 {sum += $5} END {printf "%.2f\n", sum * 100}')
    prop_gcov_4=$(cat "$base_coverage_file" | grep "genome" | awk '$2 >= 4 {sum += $5} END {printf "%.2f\n", sum * 100}')

    # Append genome coverage proportions to individual stats file
    sed -i "/^$base/s/$/\t$prop_gcov_1\t$prop_gcov_2\t$prop_gcov_3\t$prop_gcov_4/" "$individual_stats_file"

    ##############################
    # Step 9: Append individual stats to global stats file
    ##############################
    cat "$individual_stats_file" | tail -n 1 >> "$GLOBAL_STATS_OUTPUT"
done



```



# Skimming processing pipeline 1
```bash
DATAINPUT="/projects/diversea/data/20240809_Glenn_sand_eels_data/HN00220216_hdd1/RawData/"
DIR="/projects/diversea/people/sjr729/Sandeels/"
SKIM_PIP_CONSULT_KRAKEN="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/skims_processing_pipeline.sh"
conda activate skimming_echarvel
cd $DIR
DATE=$(date +%d.%m)
#module load parallel kraken2 respect consult-ii

#module load parallel kraken2 respect consult-ii
bash $SKIM_PIP_CONSULT_KRAKEN -x $DATAINPUT -r 30 -f 30 > $DIR/sandeel_skimpreprocess_$DATE.log

```
