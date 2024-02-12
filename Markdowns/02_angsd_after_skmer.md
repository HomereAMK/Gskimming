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
#!/bin/bash
#SBATCH --job-name=__BASE__FULLANGSDPROCESS
#SBATCH --output=98_log_files/__BASE__fullangsdprocees.out
#SBATCH --error=98_log_files/__BASE__fullangsdprocees.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=25g
#SBATCH --time=15:00:00
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=homerejalves.monteiro@sund.ku.dk

cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18 jdk/1.8.0_291 picard/2.27.5 parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15 python/3.11.3 openjdk/17.0.8 gatk/3.8

DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/mapped"
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skims_processing_pipeline/consult/"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
base=__BASE__

# BWA Alignment, Sorting, and Indexing
bwa mem -t "$NCPU" -R "@RG\tID:$base\tSM:$base\tPL:Illumina" "$GENOME" "$DATAINPUT"/"$base"__merged > "$DATAOUTPUT"/"$base".sam
samtools view -bS -h -q 20 -F 4 "$DATAOUTPUT"/"$base".sam > "$DATAOUTPUT"/"$base".bam
samtools sort "$DATAOUTPUT"/"$base".bam -o "$DATAOUTPUT"/"$base".sort.minq20.bam
samtools index "$DATAOUTPUT"/"$base".sort.minq20.bam
rm "$DATAOUTPUT"/"$base".sam "$DATAOUTPUT"/"$base".bam

# MarkDuplicates and Clip Overlap
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/dedup"
picard MarkDuplicates I="$DATAINPUT"/"$base".sort.minq20.bam O="$DATAOUTPUT"/"$base".nocig.dedup.minq20.bam M="$DATAOUTPUT"/"$base".duprmmetrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT
bam clipOverlap --in "$DATAOUTPUT"/"$base".nocig.dedup.minq20.bam --out "$DATAOUTPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam --stats

# Indel Realignment
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/realigned"
samtools index "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam
java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $GENOME -I "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam -o "$DATAOUTPUT"/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R $GENOME -I "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam -targetIntervals "$DATAOUTPUT"/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals -consensusDeterminationModel USE_READS -nWayOut _minq20.nocig.realigned.bam
mv *realigned.bam realigned/
mv *realigned.bai realigned/
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

