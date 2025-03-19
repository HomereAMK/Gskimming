```txt
CUNNER:
alignments: /expanse/lustre/projects/uot138/echarvel3/cunner_reads/scripts/alignments/
assembly: /expanse/lustre/projects/uot138/echarvel3/cunner_reads/genome_assemblies/GCA_020745685.1_fTauAds1.pri.cur_genomic.fna

DROSOPHILA:
alignments:  /expanse/lustre/scratch/echarvel3/temp_project/skmer2_realdata/haploid-drosophila_data/alignments/
assembly:  /expanse/lustre/scratch/echarvel3/temp_project/skmer2_realdata/haploid-drosophila_data/ncbi_dataset/data/GCA_042606445.1/GCA_042606445.1_ASM4260644v1_genomic.fna

BEE:
alignments: /expanse/lustre/scratch/echarvel3/temp_project/skmer2_realdata/haploid-bee_data/alignments/
assembly: /expanse/lustre/scratch/echarvel3/temp_project/skmer2_realdata/haploid-bee_data/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna
```

scp edcharvel@login.expanse.sdsc.edu:/expanse/lustre/scratch/echarvel3/temp_project/skmer2_realdata/haploid-bee_data/alignments/* 

#### Bee ####
```bash
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18 jdk/1.8.0_291 picard/2.27.5 parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15 python/3.11.3 openjdk/17.0.8 gatk/3.8


DATAOUTPUT="/projects/tomg/people/sjr729/echarvel_bam/bee/realigned"
DATAINPUT="/projects/tomg/people/sjr729/echarvel_bam/bee/bam/4x/dedupe_bam_files"
GENOME="/projects/tomg/people/sjr729/echarvel_bam/bee/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
GATK_JAR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"




##### Create the output directory if it doesn't exist
mkdir -p $DATAOUTPUT
samtools faidx $GENOME
picard CreateSequenceDictionary \
    R="$GENOME"\
    O="$GENOME".dict
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