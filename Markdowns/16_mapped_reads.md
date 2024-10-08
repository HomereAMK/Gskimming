## 1) raw reads Mapping to the NO_INV_genome
# Herring
```bash 
cd /projects/mjolnir1/people/sjr729/ #change path

# Load necessary modules
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18 jdk/1.8.0_291 picard/2.27.5 parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15 python/3.11.3 openjdk/17.0.8 gatk/3.8

# Define the input and output directories, and the genome reference
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/rawfastqmodernClupea/"
DATAOUTPUT="/projects/mjolnir1/people/sjr729/Skmer_ms/Herring/NoInv_RawReadsMappedToNoInvRef"
NOINV_GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic_noInversion.fna"

# Create the output directory if it doesn't exist
mkdir -p $DATAOUTPUT

# Loop through each _1.fastq.gz file in the DATAINPUT directory
for filepath1 in "$DATAINPUT"/*_1.fastq.gz; do
    # Derive the corresponding _2.fastq.gz filepath
    filepath2="${filepath1/_1.fastq.gz/_2.fastq.gz}"
    
    # Extract the base name (prefix before _1.fastq.gz)
    filename=$(basename "$filepath1")
    base=${filename%_1.fastq.gz}

    # Check if the BAM file already exists
    if [ -f "$DATAOUTPUT/$base.noinv.sort.minq20.bam" ]; then
        echo "Skipping $base, BAM file already exists."
        continue
    fi

    echo "Processing $base"

    # BWA Alignment using paired-end reads
    bwa mem -t 20 -R "@RG\tID:$base\tSM:$base\tPL:Illumina" "$NOINV_GENOME" "$filepath1" "$filepath2" > "$DATAOUTPUT/$base.sam" 2> "$DATAOUTPUT/$base.bwa.log"

    # Convert SAM to BAM, filter, and sort
    samtools view -bS -h -q 20 -F 4 "$DATAOUTPUT/$base.sam" > "$DATAOUTPUT/$base.bam"
    samtools sort "$DATAOUTPUT/$base.bam" -o "$DATAOUTPUT/$base.noinv.sort.minq20.bam"

    # Index the sorted BAM
    samtools index "$DATAOUTPUT/$base.noinv.sort.minq20.bam"

    # Remove intermediate files
    rm "$DATAOUTPUT/$base.sam" "$DATAOUTPUT/$base.bam"
done

```


## 2) bam to fastq
```bash
module load bedtools parallel
cd /projects/mjolnir1/people/sjr729/Skmer_ms/Herring/ #change path
# Create the output directory if it doesn't exist
output_dir="/projects/mjolnir1/people/sjr729/Skmer_ms/Herring/NoInv_RawReadsbbmap_mapped_fq" 
mkdir -p $output_dir
input_dir="/projects/mjolnir1/people/sjr729/Skmer_ms/Herring/NoInv_RawReadsMappedToNoInvRef"
# Change to the directory containing the .bam files
cd $input_dir

# Loop 
for bam_file in *.noinv.sort.minq20.bam; do
    fq_file="${bam_file%.bam}.fq"
    output_path="$output_dir/$fq_file"
    start_time=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$start_time] Starting conversion of $bam_file to $output_path" >> $output_dir/conversion.log
    bedtools bamtofastq -i $bam_file -fq $output_path
    end_time=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$end_time] Finished conversion of $bam_file to $output_path" >> $output_dir/conversion.log
done

```

## 3) raw reads to bbmap reads
```bash
#!/usr/bin/env bash

# Set the directory where the .fq files are located
input_dir="/projects/mjolnir1/people/sjr729/Skmer_ms/Herring/NoInv_RawReadsMappedToNoInvRef"
output_dir="/projects/mjolnir1/people/sjr729/Skmer_ms/Herring/NoInv_RawReadsMappedToNoInvRefbbmap_Processed_fq"
dukmem=16
splitsize=20000000
bbmapdir="/projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/bbmap"  # Adjust the bbmap directory accordingly

# Loop over .fq files in the input directory
for fq_file in ${input_dir}/*.fq; do
  output_file="${output_dir}/$(basename ${fq_file/.fq.gz/_merged.fq})"
  
  # Temporary file for splitting
  split=`mktemp -t`

  # Split the input .fq file
  zcat -f ${fq_file} | split -l ${splitsize} -d - ${split}_

  for x in ${split}_*; do
    # Adapter removal
    ${bbmapdir}/bbduk.sh -Xmx${dukmem}g -Xms${dukmem}g in=$x \
      out=${x}_DUK ref=adapters,phix ktrim=r k=23 mink=11 hdist=1 tpe tbo overwrite=true

    rm $x

    # Deduplication
    ${bbmapdir}/dedupe.sh -Xmx${dukmem}g -Xms${dukmem}g in=${x}_DUK \
      out=${x}_MERGED overwrite=true

    rm ${x}_DUK
  done

  # Merge the final output files into the output directory
  cat ${split}_MERGED_* > ${output_file}

  # Clean up temporary files
  rm ${split} ${split}_MERGED_*
done


