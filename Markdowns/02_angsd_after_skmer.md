## This script takes consult Skmer preprocess files and map to genome ref, remove deduplicated read and realigned reads around indels
# For AtmoreClupea
# Needs 

# Full angsd preprocess = full_launcher.sh
```bash
#!/bin/bash
#SBATCH --job-name=__BASE__FULLANGSDPROCESS
#SBATCH --output=98_log_files/__BASE__fullangsdprocees.out
#SBATCH --error=98_log_files/__BASE__fullangsdprocees.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=25g   
#SBATCH --time=15:00:00
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=homerejalves.monteiro@sund.ku.dk

#pwd
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/
#fasta seq dictionary file ref picard
#picard CreateSequenceDictionary R= $GENOME
#-fai ref
#samtools faidx $GENOME


#module
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18  jdk/1.8.0_291  picard/2.27.5  parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15 python/3.11.3 openjdk/17.0.8 gatk/3.8 

module load gcc/8.2.0 bwa/0.7.17 samtools/1.18  jdk/1.8.0_291  picard/2.27.5  parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15 python/3.11.3 openjdk/17.0.8  gatk-framework/


# Global variables
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/mapped"
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skims_processing_pipeline/consult/"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
base=__BASE__

### BWA ###
#     Align reads
    echo "Aligning $base"
    ID=$(echo "@RG\tID:$base\tSM:$base\tPL:Illumina")

  # Align reads 1 step
    bwa mem -t "$NCPU" \
        -R "$ID" \
        "$GENOME" \
        "$DATAINPUT"/"$base"__merged > "$DATAOUTPUT"/"$base".sam


        # Create bam file
    echo "Creating bam for $base"

    samtools view -bS -h -q 20 -F 4 \
    "$DATAOUTPUT"/"$base".sam >"$DATAOUTPUT"/"$base".bam


     echo "Creating sorted bam for $base"
        samtools sort "$DATAOUTPUT"/"$base".bam -o "$DATAOUTPUT"/"$base".sort.minq20.bam
        samtools index "$DATAOUTPUT"/"$base".sort.minq20.bam
   
   # Clean up
    echo "Removing "$DATAOUTPUT"/"$base".sam"
    echo "Removing "$DATAOUTPUT"/"$base".bam"

        rm "$DATAOUTPUT"/"$base".sam
        rm "$DATAOUTPUT"/"$base".bam

wait

### Mark 'n' clip duplicates ###
# Global variables
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/dedup"
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/mapped"

#tryout with NO CIGAR on MarkDuplicates
picard MarkDuplicates \
I="$DATAINPUT"/"$base".sort.minq20.bam \
O="$DATAOUTPUT"/"$base".nocig.dedup.minq20.bam \
M="$DATAOUTPUT"/"$base".duprmmetrics.txt \
REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT

wait 

#scripts ClipOverlap with NO CIGAR on MarkDuplicates
/projects/mjolnir1/apps/conda/bamutil-1.0.15/bin/bam clipOverlap \
--in "$DATAOUTPUT"/"$base".nocig.dedup.minq20.bam \
--out "$DATAOUTPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam \
--stats

wait


### Realign around indels ###
# Global variables
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/realigned"
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/dedup"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"

#move to present working dir
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/


# Index bam files
samtools index "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam 

wait 

## Create list of potential in-dels nocig
java -jar /projects/mjolnir1/apps/conda/gatk-3.8/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $GENOME \
-I "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam  \
-o "$DATAOUTPUT"/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals 

wait

## Run the indel realigner tool nocig
java -jar /projects/mjolnir1/apps/conda/gatk-3.8/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $GENOME \
-I "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam \
-targetIntervals "$DATAOUTPUT"/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals \
--consensusDeterminationModel USE_READS --nWayOut _minq20.nocig.realigned.bam

##
mv *realigned.bam realigned/
mv *realigned.bai realigned/


```

