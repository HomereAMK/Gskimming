#!/bin/bash
#SBATCH --job-name=4x__BASE__FULLANGSDPROCESS
#SBATCH --output=98_log_files/4x_BASE__fullangsdprocess.out
#SBATCH --error=98_log_files/4x_BASE__fullangsdprocess.err
#SBATCH --cpus-per-task=12
#SBATCH --mem=25g
#SBATCH --time=24:00:00
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails

#pwd
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/
#fasta seq dictionary file ref picard
#picard CreateSequenceDictionary R= $GENOME
#-fai ref
#samtools faidx $GENOME


#module
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18  jdk/1.8.0_291  picard/2.27.5  
module load parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15  python/3.11.3 openjdk/17.0.8 



# Global variables
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/4x_mapped"
#DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/subsampled_4/"
#GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
base=__BASE__

### BWA ###
#     Align reads
#    echo "Aligning $base"
#    ID=$(echo "@RG\tID:$base\tSM:$base\tPL:Illumina")

  # Align reads 1 step
#    bwa mem -t 8 \
#        -R "$ID" \
#        "$GENOME" \
#        "$DATAINPUT"/unclassified-kra_unclassified-seq_"$base"__merged.fastq > "$DATAOUTPUT"/"$base".sam


        # Create bam file
#    echo "Creating bam for $base"
#    samtools view -bS -h -q 20 -F 4 \
#    "$DATAOUTPUT"/"$base".sam >"$DATAOUTPUT"/"$base".bam


#     echo "Creating sorted bam for $base"
#        samtools sort "$DATAOUTPUT"/"$base".bam -o "$DATAOUTPUT"/"$base".sort.minq20.bam
#        samtools index "$DATAOUTPUT"/"$base".sort.minq20.bam

   # Clean up
#    echo "Removing "$DATAOUTPUT"/"$base".sam"
#    echo "Removing "$DATAOUTPUT"/"$base".bam"

#       rm "$DATAOUTPUT"/"$base".sam
#        rm "$DATAOUTPUT"/"$base".bam

#wait

### Mark 'n' clip duplicates ###
# Global variables
# Global variables
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/4x_dedup"
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/4x_mapped"

# Create the output directory if it doesn't exist
mkdir -p $DATAOUTPUT

# Loop through each .sort.minq20.bam file in the DATAINPUT directory
for filepath in "$DATAINPUT"/*.sort.minq20.bam; do
  # Extract the base name (prefix before the first ".")
  filename=$(basename "$filepath")
  base=${filename%%.*}

  echo "Processing $base"

  # Perform the MarkDuplicates operation
  picard MarkDuplicates \
    I="$DATAINPUT"/"$base".sort.minq20.bam \
    O="$DATAOUTPUT"/"$base".nocig.dedup.minq20.bam \
    M="$DATAOUTPUT"/"$base".duprmmetrics.txt \
    REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT

done


#scripts ClipOverlap with NO CIGAR on MarkDuplicates
#/projects/mjolnir1/apps/conda/bamutil-1.0.15/bin/bam clipOverlap \
#--in "$DATAOUTPUT"/"$base".nocig.dedup.minq20.bam \
#--out "$DATAOUTPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam \
#--stats

#wait


### Realign around indels ###
# Global variables
#DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/4x_realigned"
#DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/4x_dedup"
#GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"

#move to present working dir
#cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/


# Index bam files
#samtools index "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam

#wait

## Create list of potential in-dels nocig
#java -jar /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $GENOME \
#-I "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam  \
#-o "$DATAOUTPUT"/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals

#wait

## Run the indel realigner tool nocig
#java -jar /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
#-T IndelRealigner \
#-R $GENOME \
#-I "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam \
#-targetIntervals "$DATAOUTPUT"/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals \
#--consensusDeterminationModel USE_READS --nWayOut _minq20.nocig.realigned.bam

##
#mv *realigned.bam realigned/
#mv *realigned.bai realigned/


#!/bin/bash

#module
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18  jdk/1.8.0_291  picard/2.27.5  
module load parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15  python/3.11.3 openjdk/17.0.8 

# Global variables
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/4x_mapped"
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/realigned"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
GATK_JAR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar" 

# Create the output directory if it doesn't exist
mkdir -p $DATAOUTPUT

# Loop through each .nocig.dedup_clipoverlap.minq20.bam file in the DATAINPUT directory
for filepath in "$DATAINPUT"/*.nocig.dedup_clipoverlap.minq20.bam; do
  # Extract the base name (prefix before the first ".")
  filename=$(basename "$filepath")
  base=${filename%%.*}

  echo "Processing $base"

  # Index the BAM file
  samtools index "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam

  # Create realignment targets
  java -jar $GATK_JAR -T RealignerTargetCreator \
    -R $GENOME \
    -I "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam \
    -o "$DATAOUTPUT"/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals

  # Perform the realignment
  java -jar $GATK_JAR -T IndelRealigner \
    -R $GENOME \
    -I "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam \
    -targetIntervals "$DATAOUTPUT"/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals \
    -consensusDeterminationModel USE_READS \
    -nWayOut "$DATAOUTPUT"/"$base"_minq20.nocig.realigned.bam

  # Move the realigned files to the output directory
  mv "$DATAOUTPUT"/*realigned.bam "$DATAOUTPUT"
  mv "$DATAOUTPUT"/*realigned.bai "$DATAOUTPUT"

done
