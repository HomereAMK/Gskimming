
# Set up the file structure.
```bash

BASEDIR=/projects/mjolnir1/people/sjr729/Herring_4x_loco
mkdir docs
mkdir config

#Prepare a sample table 
Rscript/loco-pipe.R #change $bam path in the table


#Prepare a sample table and a chromosome table
# Chromosome table for herring dataset
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18
REFERENCE_GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
OUTPUT_FILE="/projects/mjolnir1/people/sjr729/Herring_4x_loco/docs/chr_table.tsv"
ANNOTATION_FILE="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna.ann"
# Extract the required information and create the chromosome table
awk '/chromosome:/ {split($0, arr, ":"); split(arr[2], brr, ","); print $2 "\t" brr[1]}' "$ANNOTATION_FILE" > "$OUTPUT_FILE"


