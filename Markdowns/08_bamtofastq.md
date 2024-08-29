## raw reads Mapping to the Ref genome

```bash
cd /projects/mjolnir1/people/sjr729/ #change path
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18 jdk/1.8.0_291 picard/2.27.5 parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15 python/3.11.3 openjdk/17.0.8 gatk/3.8
DATAOUTPUT="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/NoInv_bbmapReadsMappedToNoInvRef"
DATAINPUT="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/skmer1/echarvel_oedulis_skimpip_results_17.08_nonrandlib_tresh0.2_0.2/bbmap_reads/"
NOINV_GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/ncbi_dataset/data/GCF_947568905.1/GCF_947568905.1_xbOstEdul1.1_genomic_noInversion.fna"

# Create the output directory if it doesn't exist
mkdir -p $DATAOUTPUT

# Loop through each _1.fastq.gz file in the DATAINPUT directory
#for filepath1 in "$DATAINPUT"/*_read1.fq.gz; do
#  # Derive the corresponding _2.fastq.gz filepath
#  filepath2="${filepath1/_read1.fq.gz/_read2.fq.gz}"
  
  # Extract the base name (prefix before _1.fastq.gz)
#  filename=$(basename "$filepath1")
#  base=${filename%_read1.fq.gz}

  # Check if the BAM file already exists
#  if [ -f "$DATAOUTPUT/$base.sort.minq20.bam" ]; then
#    echo "Skipping $base, BAM file already exists."
#    continue
#  fi

# Loop through each fastq file in the DATAINPUT directory
for filepath in "$DATAINPUT"/*.fastq; do
  # Extract the base name (prefix after "unclassified-kra_" and before ".fq")
  filename=$(basename "$filepath")
  base=${filename}
  base=${base%%.fastq}

  # Check if the BAM file already exists
  if [ -f "$DATAOUTPUT/$base.noinv.sort.minq20.bam" ]; then
    echo "Skipping $base, BAM file already exists."
    continue
  fi

  echo "Processing $base"

  # BWA Alignment
  bwa mem -t 20 -R "@RG\tID:$base\tSM:$base\tPL:Illumina" "$NOINV_GENOME" "$filepath" > "$DATAOUTPUT/$base.sam" 2> "$DATAOUTPUT/$base.bwa.log"

  # BWA Alignment using paired-end reads
  #bwa mem -t 18 -R "@RG\tID:$base\tSM:$base\tPL:Illumina" "$GENOME" "$filepath1" "$filepath2" > "$DATAOUTPUT/$base.sam" 2> "$DATAOUTPUT/$base.bwa.log"

  # Convert SAM to BAM, filter, and sort
  samtools view -bS -h -q 20 -F 4 "$DATAOUTPUT/$base.sam" > "$DATAOUTPUT/$base.bam"
  samtools sort "$DATAOUTPUT/$base.bam" -o "$DATAOUTPUT/$base.noinv.sort.minq20.bam"

  # Index the sorted BAM
  samtools index "$DATAOUTPUT/$base.noinv.sort.minq20.bam"

  # Remove intermediate files
  rm "$DATAOUTPUT/$base.sam" "$DATAOUTPUT/$base.bam"

done


```
## Converting BAM to FASTQ with Bedtools and GNU Parallel


# Set locale settings
export LANGUAGE=en_DK.UTF-8
export LC_ALL=en_DK.UTF-8
export LC_CTYPE=en_DK.UTF-8
export LANG=en_DK.UTF-8

# module Bedtools and  Parallel
module load bedtools parallel

# code
```bash
cd /projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/ #change pathe
# Create the output directory if it doesn't exist
output_dir="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/bbmap_mapped_fq" 
mkdir -p $output_dir

# Change to the directory containing the .bam files
cd bbmapReadsMappedToRef


# Loop 
for bam_file in *.sort.minq20.bam; do
    fq_file="${bam_file%.bam}.fq"
    output_path="$output_dir/$fq_file"
    start_time=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$start_time] Starting conversion of $bam_file to $output_path" >> $output_dir/conversion.log
    bedtools bamtofastq -i $bam_file -fq $output_path
    end_time=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$end_time] Finished conversion of $bam_file to $output_path" >> $output_dir/conversion.log
done

export -f convert_bam_to_fastq

# Find all BAM files and process them in parallel
find . -name 'unclassified-seq_*.bam' | parallel convert_bam_to_fastq
```

##   wrangling to have genomes .fa without inversions
```bash
#Herring
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
CHR6="LR535862.1"
CHR12="LR535868.1"
CHR17="LR535873.1"
CHR23="LR535879.1"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"

awk '
/^>LR535862.1/ {remove=1}
/^>LR535868.1/ {remove=1}
/^>LR535873.1/ {remove=1}
/^>LR535879.1/ {remove=1}
$0 ~ /^>/ && !/^>LR535862.1/ && !/^>LR535868.1/ && !/^>LR535873.1/ && !/^>LR535879.1/ {remove=0}
remove == 0 {print}
' $GENOME | grep -v "^--$" > GCA_900700415.2_Ch_v2.0.2_genomic_noInversion.fna

#Oedulis
#chromosome 4, 5 and 8
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/ncbi_dataset/data/GCF_947568905.1/GCF_947568905.1_xbOstEdul1.1_genomic.fna"
awk '
/^>NC_079167.1/ {remove=1}
/^>NC_079168.1/ {remove=1}
/^>NC_079171.1/ {remove=1}
$0 ~ /^>/ && !/^>NC_079167.1/ && !/^>NC_079168.1/ && !/^>NC_079171.1/ {remove=0}
remove == 0 {print}
' $GENOME | grep -v "^--$" > GCF_947568905.1_xbOstEdul1.1_genomic_noInversion.fna


```

## Use run_skmer.sh to get distance matrix of skimmed preprocessed read mapped to the ref genome
```bash
conda activate skimming_echarvel
cd /projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/bbmapReadsMappedToRef/
module load parallel
DATE=$(date +%d.%m)
bash /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/run_skmer.sh -i /projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/bbmapReadsMappedToRef/ -t 19 -p 10 -o /projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/skmer1/library_runskmer_bbmapreads_mappedtoRef_skmer1_Oedulis_$DATE
 # result are in library_runskmer_mappedreads_skmer1_condaskmer2test_26.06
 ```

## Skmer2 Distance Calculation with skmer2 

```bash
#Herring bbmap reads map to ref 
cd /projects/mjolnir1/people/sjr729/Skmer_ms/Herring/bbmapReadsMappedToRef/
conda activate skimming_echarvel
#conda install jellyfish seqtk mash
module load parallel
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
DATE=$(date +%d.%m)
# build the skmer2 library
python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/__main_TESTING.py --debug reference /projects/mjolnir1/people/sjr729/Skmer_ms/Herring/bbmap_mapped_fq/ -r $GENOME -p 2 -o "library_skmer2_bbmapreads_mappedtoRef_Herring_$DATE-dist-mat" -l library_skmer2_bbmapreads_mappedtoRef_Herring_$DATE
