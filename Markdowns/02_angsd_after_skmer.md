## Skmer Preprocessing Mapping and Realignment for AtmoreClupea

This script processes Skmer preprocess files, maps them to a reference genome, removes deduplicated reads, and realigns reads around indels.

### Create a List of Merged Reads for Clupea Atmore with Full Path

```bash
DATE=$(date +%Y-%m-%d)  # Sets DATE to the format YYYY-MM-DD
ls -d unclassified-kra_* > "${DATE}_list_Atmore_mergedreads.txt"
```

### Extract IDs and Construct Filename Patterns

```bash
while IFS= read -r line; do
    id_part=$(echo "$line" | grep -o 'ERR97092[0-9]*')
    find "/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skims_processing_pipeline/consult/" -type f -name "*${id_part}*" >> Atmore_modern_full_paths.txt
done < "2024-02-03_list_Atmore_mergedreads.txt"
```

### Move Non-matching Files to a New Directory

```bash
#!/bin/bash
SOURCE_DIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skims_processing_pipeline/consult/"
TARGET_DIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skims_processing_pipeline/consult_nonmodern/"
LIST_FILE="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/Atmore_modern_full_paths.txt"
mkdir -p "$TARGET_DIR"
find "$SOURCE_DIR" -type f | while read -r file; do
    if ! grep -Fxq "$file" "$LIST_FILE"; then
        mv "$file" "$TARGET_DIR"
    fi
done
```


### Full ANGSD Preprocessing: Full Launcher Script

```bash
### mapping 
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18 jdk/1.8.0_291 picard/2.27.5 parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15 python/3.11.3 openjdk/17.0.8 gatk/3.8

DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/1x_mapped"
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/subsampled_1/"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"

# Create the output directory if it doesn't exist
mkdir -p $DATAOUTPUT

# Loop through each merged fastq file in the DATAINPUT directory
for filepath in "$DATAINPUT"/*__merged.fastq; do
  # Extract the base name (prefix after "unclassified-kra_unclassified-seq_" and before "__merged")
  filename=$(basename "$filepath")
  base=${filename#unclassified-kra_unclassified-seq_}
  base=${base%%__merged.fastq}

  # Check if the BAM file already exists
  if [ -f "$DATAOUTPUT/$base.sort.minq20.bam" ]; then
    echo "Skipping $base, BAM file already exists."
    continue
  fi

  echo "Processing $base"

  # BWA Alignment
  bwa mem -t 8 -R "@RG\tID:$base\tSM:$base\tPL:Illumina" "$GENOME" "$DATAINPUT"/"$filename" > "$DATAOUTPUT"/"$base".sam 2> "$DATAOUTPUT"/"$base".bwa.log

  # Convert SAM to BAM, filter, and sort
  samtools view -bS -h -q 20 -F 4 "$DATAOUTPUT"/"$base".sam > "$DATAOUTPUT"/"$base".bam
  samtools sort "$DATAOUTPUT"/"$base".bam -o "$DATAOUTPUT"/"$base".sort.minq20.bam

  # Index the sorted BAM
  samtools index "$DATAOUTPUT"/"$base".sort.minq20.bam

  # Remove intermediate files
  rm "$DATAOUTPUT"/"$base".sam "$DATAOUTPUT"/"$base".bam

done

### MarkDuplicates and Clip Overlap
# Global variables
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/2x_dedup" #change path
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/2x_mapped" #change path

# Create the output directory if it doesn't exist
mkdir -p $DATAOUTPUT

# Loop through each .sort.minq20.bam file in the DATAINPUT directory
for filepath in "$DATAINPUT"/*.sort.minq20.bam; do
  # Extract the base name (prefix before the first ".")
  filename=$(basename "$filepath")
  base=${filename%%.*}

  # Check if the deduplicated BAM file already exists
  if [ -f "$DATAOUTPUT/$base.nocig.dedup.minq20.bam" ]; then
    echo "Skipping $base, deduplicated BAM file already exists."
    continue
  fi

  echo "Processing $base"

  # Perform the MarkDuplicates operation
  picard MarkDuplicates \
    I="$DATAINPUT"/"$base".sort.minq20.bam \
    O="$DATAOUTPUT"/"$base".nocig.dedup.minq20.bam \
    M="$DATAOUTPUT"/"$base".duprmmetrics.txt \
    REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT

done




# Indel Realignment
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18 jdk/1.8.0_291 picard/2.27.5
module load parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15 python/3.11.3 openjdk/17.0.8 gatk/3.8

# Global variables
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/2x_dedup" #change path
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/2x_realigned" #change path
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
GATK_JAR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"

# Create the output directory if it doesn't exist
mkdir -p $DATAOUTPUT

# Loop through each .bam file in the DATAINPUT directory
for filepath in "$DATAINPUT"/*.nocig.dedup.minq20.bam; do
  # Extract the base name (prefix before the first ".")
  filename=$(basename "$filepath")
  base=${filename%%.*}

  # Check if the realigned BAM file already exists
  if [ -f "$DATAOUTPUT/$base"_minq20.nocig.realigned.bam ]; then
    echo "Skipping $base, realigned BAM file already exists."
    continue
  fi

  echo "Processing $base"

  # Index the BAM file
  samtools index "$DATAINPUT"/"$base".nocig.dedup.minq20.bam

  # Create realignment targets
  java -jar $GATK_JAR -T RealignerTargetCreator \
    -R $GENOME \
    -I "$DATAINPUT"/"$base".nocig.dedup.minq20.bam \
    -o "$DATAOUTPUT"/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals

  # Perform the realignment
  java -jar $GATK_JAR -T IndelRealigner \
    -R $GENOME \
    -I "$DATAINPUT"/"$base".nocig.dedup.minq20.bam \
    -targetIntervals "$DATAOUTPUT"/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals \
    -model USE_READS \
    -o "$DATAOUTPUT"/"$base"_minq20.nocig.realigned.bam

done

# Move the realigned files to the output directory
mv "$DATAOUTPUT"/*realigned.bam "$DATAOUTPUT"
mv "$DATAOUTPUT"/*realigned.bai "$DATAOUTPUT"
```




















### Create the BAM List for Modern Clupea Atmore

```bash
SEARCH_DIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/realigned/"
LIST_FILE="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/Atmore_modern_full_paths.txt"
OUTPUT_FILE="Atmore_modern_bam_list.txt"
BAMS=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/realigned/*.bam
ls $BAMS > 01_infofiles/Atmore_modern_bam_list.txt
```

### Downsizing bam files
```bash
awk -F, '{print "reformat.sh in1=done/"$1"_1.fastq.gz in2=done/"$1"_2.fastq.gz out1="$1"_4x_1.fastq.gz out2="$1"_4x_2.fastq.gz samplerate="$6"\njava -ea -Xms300m -cp /projects/mjolnir1/apps/conda/bbmap-39.01/opt/bbmap-39.01-0/current/ jgi.ReformatReads in1=done/"$1"_1.fastq.gz in2=done/"$1"_2.fastq.gz out1="$1"_4x_1.fastq.gz out2="$1"_4x_2.fastq.gz samplerate="$6"}' library_4x/feb24_Atmore_modern_4x_stat_clupea_df.csv
```



### Stats mapping Clupea 
module load tools
module load ngs
module load samtools/1.14

```bash
DIR=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea

# Iterate over only the *_1.fastq.gz files to avoid processing each base twice
for file in $DIR/rawfastqmodernClupea/*_1.fastq.gz; do
    # Extract the base name by taking the first 10 characters of the filename
    base=$(basename $file _1.fastq.gz | cut -c1-10)
    #Global variables
    #raw reads
    a=`zcat $DIR/rawfastqmodernClupea/"$base"_1.fastq.gz  | wc -l | awk '{print $1/4}'` #raw read forward
    b=`zcat $DIR/rawfastqmodernClupea/"$base"_2.fastq.gz | wc -l | awk '{print $1/4}'` #raw read reverse
    echo $(( $a + $b )) > $DIR/StatsMapping/"$base".count_fastq_1.tmp
    #raw bases
    c=`zcat $DIR/rawfastqmodernClupea/"$base"_1.fastq.gz | awk 'NR%4==2' | tr -d "\n" | wc -m` 
    d=`zcat $DIR/rawfastqmodernClupea/"$base"_2.fastq.gz | awk 'NR%4==2' | tr -d "\n" | wc -m`
    echo $(( $c + $d )) > $DIR/StatsMapping/"$base".count_fastq_2.tmp

    #trim bases
    e=`cat $DIR/skims_processing_pipeline_jan24/consult/unclassified-seq_"$base"__merged | awk 'NR%4==2' | tr -d "\n" | wc -m` 
    echo  $e  > $DIR/StatsMapping/"$base".count_fastq_3.tmp 

    #mapped bases
    samtools stats $DIR/angsd/mapped/unclassified-seq_"$base".sort.minq20.bam -@ 12 | grep ^SN | cut -f 2- | grep "^bases mapped (cigar)" | cut -f 2 > $DIR/StatsMapping/"$base".count_bam_1.tmp
        
    #deduplicate mapped bases
    samtools stats $DIR/angsd/dedup/unclassified-seq_"$base".nocig.dedup_clipoverlap.minq20.bam -@ 12 | grep ^SN | cut -f 2- | grep "^bases mapped (cigar)" | cut -f 2  > $DIR/StatsMapping/"$base".count_bam_2.tmp

    #realigned around indels mapped bases
    samtools stats $DIR/angsd/realigned/unclassified-seq_"$base".nocig.dedup_clipoverlap.minq20_minq20.nocig.realigned.bam -@ 12 | grep ^SN | cut -f 2- | grep "^bases mapped (cigar)" | cut -f 2  > $DIR/StatsMapping/"$base".count_bam_3.tmp

    RAWREADS=`cat $DIR/StatsMapping/"$base".count_fastq_1.tmp`
    RAWBASES=`cat $DIR/StatsMapping/"$base".count_fastq_2.tmp`
    ADPTERCLIPBASES=`cat $DIR/StatsMapping/"$base".count_fastq_3.tmp`
    MAPPEDBASES=`cat $DIR/StatsMapping/"$base".count_bam_1.tmp`
    DEDUPMAPPEDBASES=`cat $DIR/StatsMapping/"$base".count_bam_2.tmp`
    REALIGNEDMAPPEDBASES=`cat $DIR/StatsMapping/"$base".count_bam_3.tmp`

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $base $RAWREADS $RAWBASES $ADPTERCLIPBASES $MAPPEDBASES $DEDUPMAPPEDBASES $REALIGNEDMAPPEDBASES >> $DIR/StatsMapping/Summary_modernClupea_stats_9apr23.txt
done
```



Mapping to the Reference Genome:
Alignment, Sorting, and Indexing: Reads were aligned to the reference genome using BWA, then sorted and indexed with Samtools.
Duplicate Marking and Overlap Clipping: Duplicates were marked using Picard's MarkDuplicates, and overlaps were clipped using BAMUtil.
Indel Realignment: Realignment around indels was performed using GATK's RealignerTargetCreator and IndelRealigner.
Variant Site Detection:

Variants were detected in the population with a p-value threshold of ≤ 1e-6.
High Linkage Disequilibrium Removal:

Regions with high linkage disequilibrium were excluded. Specifically, SNPs on chromosomes 6, 12, 17, and 23 were excluded due to known inversion regions (Jamsandekar et al., 2023).





Mapping to the Reference Genome:
Align reads to the reference genome using BWA.
Sort and index the aligned reads with Samtools.
Mark duplicates using Picard's MarkDuplicates.
Clip overlaps with BAMUtil.
Perform indel realignment with GATK's RealignerTargetCreator and IndelRealigner.
Variant Site Detection:

Detect variants with a p-value threshold of ≤ 1e-6.
High Linkage Disequilibrium Removal: Exclude regions with high linkage disequilibrium.
Remove SNPs on chromosomes 6, 12, 17, and 23 due to known inversion regions (Jamsandekar et al., 2023).


Dataset: Atlantic herring dataset downloaded from ENA/SNCBI (Atmore et al., 2022).

Species Information: Atlantic herring is one of the most abundant vertebrates on Earth with an effective population size (Ne) over a million and a census population size (Nc) over a trillion. The species has adapted to various ecological and environmental conditions such as variations in salinity, water temperature, light conditions, spawning seasons, and food resources.

Population Structure Analysis:
The analysis was performed using an MDS plot with a distance matrix obtained from state-of-the-art software for genome skims (low-coverage whole genome sequencing data, lcWGS): ANGSD and Skmer.The dataset includes 45 modern individuals from 9 different sampling sites sequenced below 10x coverage.