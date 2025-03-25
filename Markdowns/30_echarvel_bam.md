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
#1.IndelRealigner
```bash
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18 jdk/1.8.0_291 picard/2.27.5 parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15 python/3.11.3 openjdk/17.0.8 gatk/3.8


DATAOUTPUT="/projects/tomg/people/sjr729/echarvel_bam/bee/6x/6x_realigned"
DATAINPUT="/projects/tomg/people/sjr729/echarvel_bam/bee/6x/dedupe_bam_files"
GENOME="/projects/tomg/people/sjr729/echarvel_bam/bee/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
GATK_JAR="/projects/tomg/people/sjr729/echarvel_bam/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"

#GENOME="/projects/tomg/people/sjr729/diversea/Skmer_ms/Cunner/genome_assemblies/GCA_020745685.1_fTauAds1.pri.cur_genomic.fna"
#DATAINPUT="/scratch03/echarvel/cunner_bams/4x_downsampled/"
#DATAOUTPUT="/scratch03/echarvel/cunner_bams/4x_realigned"
#GENOME="/scratch03/echarvel/ncbi_dataset/data/GCA_020745685.1/GCA_020745685.1_fTauAds1.pri.cur_genomic.fna"
#GATK_JAR="/calab_data/mirarab/home/hmon/softwares/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"





##### Create the output directory if it doesn't exist
mkdir -p $DATAOUTPUT
#samtools faidx $GENOME
#/calab_data/mirarab/home/echarvel3/programs/samtools-1.21

picard CreateSequenceDictionary \
    R="$GENOME"\
    O="$GENOME".dict
# Loop through each .bam file in the DATAINPUT directory
for filepath in "$DATAINPUT"/*.nocig.dedup.minq20.bam; do
  # Extract the base name (prefix before the first ".")
  filename=$(basename "$filepath")
  base=${filename%%.*}

  # Check if the realigned BAM file already exists
  if [ -f "$DATAOUTPUT/$base".nocig.dedup.minq20.bam ]; then
    echo "Skipping $base, realigned BAM file already exists."
    continue
  fi

  echo "Processing $base"

  # Index the BAM file
  samtools index "$DATAINPUT"/"$base".nocig.dedup.minq20.bam
  #/calab_data/mirarab/home/echarvel3/programs/samtools-1.21/samtools index "$DATAINPUT"/"$base".nocig.dedup.minq20.bam

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

#2.Loco-pipe with -isHap and without for heterozygosity

```bash
## Create the sample_table tsv
cd /projects/tomg/people/sjr729/echarvel_bam/bee

# Create the header
echo -e "sample_name\tbam\tpopulation" > samples_no_species.tsv

# Append lines for each BAM
for f in 4x_realigned/unclassified-seq_*_minq20.nocig.realigned.bam
do
    # Strip directory to keep only the filename:
    fname=$(basename "$f")

    # Remove prefix 'unclassified-seq_' and suffix '_minq20.nocig.realigned.bam'
    sname=$(echo "$fname" \
            | sed 's/^unclassified-seq_//' \
            | sed 's/_minq20\.nocig\.realigned\.bam$//')

    # Append row to TSV
    # Adjust 'unclassified' below if you want a different population name
    echo -e "${sname}\t$(realpath "$f")\tunclassified"
done >> samples_no_species.tsv


# 1. Extract and preserve the header from samples_no_species.tsv
head -n1 samples_no_species.tsv > header.txt

# 2. Remove the header from samples_no_species.tsv, then sort on the first column
tail -n +2 samples_no_species.tsv | sort -k1,1 > samples_no_species_sorted.tsv

# 3. Sort the colonies.txt file
sort colonies.txt > colonies_sorted.txt

# 4. Convert space-delimited lines in colonies_sorted.txt to tab-delimited
tr ' ' '\t' < colonies_sorted.txt > colonies_sorted_tab.txt

# 5. Join the two sorted files on the first column, using tabs
join -t $'\t' -1 1 -2 1 \
    samples_no_species_sorted.tsv \
    colonies_sorted_tab.txt \
    > joined_data.tsv

# 6. Add a new header line (with "colony") and append the joined data
echo -e "sample_name\tbam\tpopulation\tcolony" > Bee_colonies_table.tsv
cat joined_data.tsv >> Bee_colonies_table.tsv
```

```bash
#Run the snakemake pipeline
# Snakemake run of the loco_pipe
```bash
#Launching the pipeline
#for dry run add -n

conda activate loco-pipe
module purge
SOFTWARE_DIR=/projects/tomg/people/sjr729/diversea/
BASEDIR=/projects/tomg/people/sjr729/echarvel_bam/bee/4x_locopipe
cd $SOFTWARE_DIR
snakemake \
--use-conda \
--conda-frontend mamba \
--directory $BASEDIR \
--rerun-triggers mtime \
--scheduler greedy \
--printshellcmds \
--snakefile $SOFTWARE_DIR/loco-pipe/workflow/pipelines/loco-pipe.smk \
--cores 10 --rerun-incomplete --keep-going 

```


### CUNNER
