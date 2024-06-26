
# Set up the file structure.
```bash
SOFTWARE_DIR=/projects/mjolnir1/people/sjr729/
BASEDIR=/projects/mjolnir1/people/sjr729/Herring_raw_loco
mkdir $BASEDIR
cd $BASEDIR
mkdir docs
mkdir config
conda activate loco-pipe
#Prepare a sample table 
Rscript/loco-pipe.R #change $bam path in the table
scp with filezilla 01_infofiles/sample_table_Herring_raw_loco.tsv

#Prepare a sample table and a chromosome table
# Chromosome table for herring dataset scp /projects/mjolnir1/people/sjr729/Herring_4x_loco/docs/chr_table.tsv ./docs
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18
REFERENCE_GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
OUTPUT_FILE="/projects/mjolnir1/people/sjr729/Herring_2x_loco/docs/chr_table.tsv"
ANNOTATION_FILE="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna.ann"
# Extract the required information and create the chromosome table
awk '/chromosome:/ {split($0, arr, ":"); split(arr[2], brr, ","); print $2 "\t" brr[1]}' "$ANNOTATION_FILE" > "$OUTPUT_FILE"

#Prepare the config.yaml
/loco-pipe/config.yaml
scp with filezilla

#Plot the pipeline flowchart
BASEDIR=/projects/mjolnir1/people/sjr729/Herring_2x_loco
mkdir -p $SOFTWARE_DIR/loco-pipe/toyfish/figures/flowchart
snakemake -n --forceall --rulegraph \
--directory $BASEDIR \
--snakefile $SOFTWARE_DIR/loco-pipe/workflow/pipelines/loco-pipe.smk  | \
dot -Tpdf > $SOFTWARE_DIR/loco-pipe/toyfish/figures/flowchart/toyfish.pdf


#Launching the pipeline
#for dry run add -n
conda activate loco-pipe
SOFTWARE_DIR=/projects/mjolnir1/people/sjr729/
BASEDIR=/projects/mjolnir1/people/sjr729/Herring_raw_loco
module load mamba 
module load gcc/13.2.0
module load openjdk/20.0.0
module load R
#module load snakemake
snakemake \
--use-conda \
--conda-frontend mamba \
--directory $BASEDIR \
--rerun-triggers mtime \
--scheduler greedy \
--snakefile $SOFTWARE_DIR/loco-pipe/workflow/pipelines/loco-pipe.smk \
--cores 15
