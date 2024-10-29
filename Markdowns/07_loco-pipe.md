# Installation
```bash


SOFTWARE_DIR=/projects/mjolnir1/people/sjr729/
mamba env create -f $SOFTWARE_DIR/loco-pipe/workflow/envs/loco-pipe2.yaml 

#in lostruct env
conda activate lostruct_lcpipe
R
withr::with_libpaths(new = "/home/sjr729/miniforge3/envs/lostruct_lcpipe/lib/R/library",
                     devtools::install_github("petrelharp/local_pca/lostruct"))
```
# Set up the file structure.
```bash

SOFTWARE_DIR=/projects/mjolnir1/people/sjr729/
BASEDIR=/projects/mjolnir1/people/sjr729/EUostrea_molecol_loco
mkdir $BASEDIR
cd $BASEDIR
mkdir docs
mkdir config
conda activate loco-pipe
#Prepare a sample table 
Rscript/loco-pipe.R #change $bam path in the table
scp with filezilla 01_infofiles/SteinOedulis_loco.tsv

#Prepare a sample table and a chromosome table
# Chromosome table for herring dataset scp /projects/mjolnir1/people/sjr729/Herring_4x_loco/docs/chr_table.tsv ./docs
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18
#REFERENCE_GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
REFERENCE_GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/ncbi_dataset/data/GCF_947568905.1/GCF_947568905.1_xbOstEdul1.1_genomic.fna"
OUTPUT_FILE="/projects/mjolnir1/people/sjr729/SteinOedulis_loco/docs/chr_table.tsv"
ANNOTATION_FILE="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/ncbi_dataset/data/GCF_947568905.1/GCF_947568905.1_xbOstEdul1.1_genomic.fna.ann"
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
```
# Snakemake run of the loco_pipe
```bash
#Launching the pipeline
#for dry run add -n
conda activate loco-pipe
module purge
SOFTWARE_DIR=/projects/mjolnir1/people/sjr729/
BASEDIR=/projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2
cd $SOFTWARE_DIR
snakemake \
--use-conda \
--conda-frontend mamba \
--directory $BASEDIR \
--rerun-triggers mtime \
--scheduler greedy \
--printshellcmds \
--snakefile $SOFTWARE_DIR/loco-pipe/workflow/pipelines/loco-pipe.smk \
--cores 19 --rerun-incomplete --keep-going --unlock
#--conda-prefix /home/sjr729/miniforge3/envs \
```

# Rerun Angsd on combined SNP list to get combine Ibs.mat oedulis dataset
```bash
conda activate angsd_lcpipe
module purge

SNPLIST=/projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2/angsd/snp_calling_global/combined.subsetted.snp_list
REF=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/ncbi_dataset/data/GCF_947568905.1/GCF_947568905.1_xbOstEdul1.1_genomic.fna

angsd -bam /projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2/docs/bamlist.txt \
-ref $REF \
-sites $SNPLIST \
-P 20 \
-out /projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2/angsd/snp_calling_global/25.07_SteinOedulis_loco_run2_totalSNPlist_forDistMat \
-GL 1 \
-doGlf 2 \
-doIBS 1 \
-makematrix 1 \
-doMajorMinor 3 \
-doDepth 1 -doCounts 1 -maxDepth 100000 -dumpCounts 1 -setMinDepth 69 -setMaxDepth 230 -minInd 46 -setMinDepthInd 1 \
-minQ 20 -minMapQ 20 \
-remove_bads 1 \
-only_proper_pairs 1
```
