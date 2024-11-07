# 1. ERDA

INPUT_DIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/bam"


# 2. make var for loco-pipe

BASEDIR=/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/locopipe_hetonly
mkdir -p $BASEDIR
cd $BASEDIR
mkdir docs
mkdir config
conda activate loco-pipe
#Prepare a sample table, full path bam $1st column, pop 2nd column
INPUT_DIR_BAM="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/bam"
TANGUY_GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/Tanguy/fileOegenome10scaffoldC3G.fasta"


REFERENCE_GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/Tanguy/fileOegenome10scaffoldC3G.fasta"
OUTPUT_FILE="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/locopipe/docs/chr_table.tsv"
ANNOTATION_FILE="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/Tanguy/fileOegenome10scaffoldC3G.fasta.ann"
# manually did the chr_table_ATgenome.tsv

bam_dir="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/bam"
output_file="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/locopipe/docs/sample_table_ATgenome.tsv"
echo -e "sample_name\tspecies\tpopulation\tbam" > $output_file
for bam_file in "$bam_dir"/*.bam; do
    filename=$(basename "$bam_file")
    sample_name=${filename:0:7}
    population=${filename:0:4}
    bam_full_path="$bam_file"

    if [ "$population" == "LURI" ]; then
        species="Ostrea_lurida"
    else
        species="Ostrea_edulis"
    fi

    echo -e "${sample_name}\t${species}\t${population}\t${bam_full_path}" >> $output_file
done


#Prepare the config.yaml
/loco-pipe/config.yaml
scp with filezilla


BASEDIR=/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/locopipe_hetonly
# Snakemake run of the loco_pipe
```bash
#Launching the pipeline
#for dry run add -n
conda activate loco-pipe
module purge
SOFTWARE_DIR=/projects/mjolnir1/people/sjr729/
BASEDIR=/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/locopipe_hetonly
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


