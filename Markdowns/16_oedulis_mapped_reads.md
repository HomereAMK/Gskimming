#!/bin/bash

# Load necessary modules
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18 jdk/1.8.0_291 picard/2.27.5 parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15 python/3.11.3 openjdk/17.0.8 gatk/3.8 bedtools bbmap

# Define directories and genome reference
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/oedulis_RawFastq_14.08.24"
DATAOUTPUT="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/NoInv_RawReadsMappedToNoInvRef"
NOINV_GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/ncbi_dataset/data/GCF_947568905.1/GCF_947568905.1_xbOstEdul1.1_genomic_noInversion.fna"
FASTQ_OUTPUT="/projects/mjolnir1/people/sjr729/Skmer_ms/Herring/NoInv_RawReadsbbmap_mapped_fq"
BBMAP_OUTPUT="/projects/mjolnir1/people/sjr729/Skmer_ms/Herring/NoInv_RawReadsMappedToNoInvRefbbmap_Processed_fq"
BBMAP_DIR="/projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/bbmap"

# Create necessary output directories
mkdir -p "$DATAOUTPUT" "$FASTQ_OUTPUT" "$BBMAP_OUTPUT"

# Memory and split parameters for bbmap
dukmem=16
splitsize=20000000

# Loop through each _1.fastq.gz file in the DATAINPUT directory
for filepath1 in "$DATAINPUT"/*_1.fastq.gz; do
    filepath2="${filepath1/_1.fastq.gz/_2.fastq.gz}"
    filename=$(basename "$filepath1")
    base=${filename%_1.fastq.gz}

    # Check if BAM file already exists and skip if it does
    if [ -f "$DATAOUTPUT/$base.noinv.sort.minq20.bam" ]; then
        echo "Skipping $base, BAM file already exists."
        continue
    fi

    ##############################
    # Step 1: Raw Reads Mapping to the NO_INV_genome
    ##############################
    echo "Step 1: Mapping raw reads to NO_INV_genome for $base"
    bwa mem -t 30 -R "@RG\tID:$base\tSM:$base\tPL:Illumina" "$NOINV_GENOME" "$filepath1" "$filepath2" > "$DATAOUTPUT/$base.sam" 2> "$DATAOUTPUT/$base.bwa.log"
    
    if [ $? -ne 0 ]; then
        echo "BWA failed for $base"
        exit 1
    fi

    samtools view -bS -h -q 20 -F 4 "$DATAOUTPUT/$base.sam" > "$DATAOUTPUT/$base.bam"
    samtools sort "$DATAOUTPUT/$base.bam" -o "$DATAOUTPUT/$base.noinv.sort.minq20.bam"
    samtools index "$DATAOUTPUT/$base.noinv.sort.minq20.bam"
    rm "$DATAOUTPUT/$base.sam" "$DATAOUTPUT/$base.bam"

    if [ $? -ne 0 ]; then
        echo "SAM/BAM processing failed for $base"
        exit 1
    fi

    ##############################
    # Step 2: BAM to Fastq conversion
    ##############################
    echo "Step 2: Converting BAM to Fastq for $base"
    fq_file="$FASTQ_OUTPUT/$base.fq"

    bedtools bamtofastq -i "$DATAOUTPUT/$base.noinv.sort.minq20.bam" -fq "$fq_file"

    if [ $? -ne 0 ]; then
        echo "Conversion from BAM to Fastq failed for $base"
        exit 1
    fi

    ##############################
    # Step 3: Processing raw reads with bbmap
    ##############################
    echo "Step 3: Processing Fastq with bbmap for $base"
    split=$(mktemp -t)
    output_file="${BBMAP_OUTPUT}/${base}_merged.fq"

    zcat -f "$fq_file" | split -l "$splitsize" -d - "${split}_"

    for x in "${split}_*"; do
        # Adapter removal with bbduk
        ${BBMAP_DIR}/bbduk.sh -Xmx${dukmem}g -Xms${dukmem}g in="$x" out="${x}_DUK" ref=adapters,phix ktrim=r k=23 mink=11 hdist=1 tpe tbo overwrite=true
        if [ $? -ne 0 ]; then
            echo "bbduk failed for $x"
            exit 1
        fi

        rm "$x"

        # Deduplication with dedupe
        ${BBMAP_DIR}/dedupe.sh -Xmx${dukmem}g -Xms${dukmem}g in="${x}_DUK" out="${x}_MERGED" overwrite=true
        if [ $? -ne 0 ]; then
            echo "dedupe failed for $x"
            exit 1
        fi

        rm "${x}_DUK"
    done

    # Merge split files into final output file
    cat "${split}_MERGED_*" > "$output_file"
    rm "${split}" "${split}_MERGED_*"

    if [ $? -ne 0 ]; then
        echo "Merging failed for $fq_file"
        exit 1
    fi

    echo "Processing for $base completed successfully."
done

echo "All samples processed successfully."
