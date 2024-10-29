# 1. Rsync data (not used)
```bash
#!/bin/bash

# Define variables
REMOTE_USER="echarvel3"
REMOTE_HOST="login.expanse.sdsc.edu"
REMOTE_BASE="/expanse/lustre/projects/uot138/echarvel3/cunner_reads"
LOCAL_BASE="/projects/mjolnir1/people/sjr729/Skmer_ms/Cunner"
POP_FILES_DIR="$LOCAL_BASE/candidates_for_angsd"

# Loop over each population file
for popfile in "$POP_FILES_DIR"/*_above10x.txt; do
    # Extract the population code from the filename
    POP=$(basename "$popfile" | cut -c1-4)
    echo "Processing population: $POP"

    # Read the sample IDs from the file
    samples=($(awk '{print $1}' "$popfile"))
    
    # Randomly select 15 samples
    selected_samples=($(shuf -e "${samples[@]}" -n 15))
    
    # Loop over each selected sample
    for sample in "${selected_samples[@]}"; do
        # Construct the source and destination paths
        src="$REMOTE_USER@$REMOTE_HOST:$REMOTE_BASE/decontam/${sample}.fastq"
        dest="$LOCAL_BASE/${POP}_${sample}.fastq"
        
        # Rsync the file with the population prefix added to the filename
        rsync -av "$src" "$dest"
    done
done

```
# 2. NGDIST on published BEAGLE FILE
```bash
#Cunner
module load ngsdist
BASE_DIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Cunner"
BEAGLE="$BASE_DIR/cunner_merged.beagle.gz"
LABEL="/projects/mjolnir1/people/sjr729/Skmer_ms/Cunner/cunner_labels.tsv"
DOWNSAMP_BEAGLE="$BASE_DIR/cunner_merged_subsampled_1out50.beagle.gz"
#Create the label files on local
zcat Downloads/cunner_merged.beagle.gz | head -1 | grep -oE 'Ind[0-9]+' | sort -u > cunner_labels.tsv
#on mj 
zcat $BEAGLE | head -1 | grep -oE 'Ind[0-9]+' | awk '!seen[$0]++' > cunner_labels.tsv
#number of sites
zcat $BEAGLE | tail -n +2 | wc -l
zcat $DOWNSAMP_BEAGLE | tail -n +2 | wc -l
#Prepare a geno file by subsampling one SNP in every 50 SNPs in the beagle file
zcat $BEAGLE | awk 'NR % 50 == 0' | cut -f 4- | gzip  > $BASE_DIR/cunner_merged_subsampled_1out50.beagle.gz

#ngsdist
ngsDist --n_threads 30 --geno $DOWNSAMP_BEAGLE --pairwise_del --seed 3 --probs --n_ind 850 --n_sites 231488  --labels $LABEL --out $BASE_DIR/cunner_merged_subsampled_1out50.dist

#Performs MDS using get_PCA.R:
conda activate r_lcpipe
tail -n +3 $BASE_DIR/cunner_merged_subsampled_1out50.dist | Rscript --vanilla --slave /projects/mjolnir1/people/sjr729/get_PCA.R --no_header --data_symm -n 10 -m mds -o $BASE_DIR/cunner_merged_subsampled_1out50.mds
# ngsdist with --avg_nuc_dist
ngsDist --n_threads 30 --geno $DOWNSAMP_BEAGLE --pairwise_del --avg_nuc_dist --seed 3 --probs --n_ind 850 --n_sites 231488  --labels $LABEL --out $BASE_DIR/cunner_merged_subsampled_1out50_avg_nuc_dist.dist
tail -n +3 $BASE_DIR/cunner_merged_subsampled_1out50_avg_nuc_dist.dist | Rscript --vanilla --slave /projects/mjolnir1/people/sjr729/get_PCA.R --no_header --data_symm -n 10 -m mds -o $BASE_DIR/cunner_merged_subsampled_1out50_avg_nuc_dist.mds


#for oedulis
BASE_DIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis"
BEAGLE="$BASE_DIR/30jan23_prunedLDminweight0.5_PopStruct.beagle_noheader.gz"
LABEL="$BASE_DIR/bamlist_EUostrea.labels"
ngsDist --n_threads 30 --geno $BEAGLE --pairwise_del --seed 3 --probs --n_ind 581 --n_sites 1404180  --labels $LABEL --out $BASE_DIR/Oedulis_beagle.dist

ngsDist --n_threads 30 --geno $BEAGLE --pairwise_del --seed 3 --avg_nuc_dist --probs --n_ind 581 --n_sites 1404180  --labels $LABEL --out $BASE_DIR/Oedulis_beagle_avg_nuc_dist.dist

#mds
conda activate r_lcpipe
tail -n +3 $BASE_DIR/Oedulis_beagle.dist | Rscript --vanilla --slave /projects/mjolnir1/people/sjr729/get_PCA.R --no_header --data_symm -n 10 -m mds -o $BASE_DIR/Oedulis_beagle.mds

```


# 3. Verify Cunner BEAGLE FILE 
```bash
C_BEAGLE="/projects/mjolnir1/people/sjr729/Skmer_ms/Cunner/cunner_merged_subsampled_1out50.beagle.gz"
O_BEAGLE="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/30jan23_prunedLDminweight0.5_PopStruct.beagle.gz"

#check same format of these beagle files
zcat "/projects/mjolnir1/people/sjr729/Skmer_ms/Cunner/cunner_merged_subsampled_1out50.beagle.gz" | head -n 1 # this file got an header
zcat "/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/30jan23_prunedLDminweight0.5_PopStruct.beagle.gz" | head -n 1 #no header

#do a noheader version of C_BEAGLE
C_BEAGLE="/projects/mjolnir1/people/sjr729/Skmer_ms/Cunner/cunner_merged_subsampled_1out50.beagle.gz"
zcat $C_BEAGLE | tail -n +2 > /projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/30jan23_prunedLDminweight0.5_PopStruct.beagle_noheader.gz
BEAGLE_NOHEAD="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/30jan23_prunedLDminweight0.5_PopStruct.beagle_noheader.gz"
BASE_DIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Cunner"
ngsDist --n_threads 30 --geno $BEAGLE_NOHEAD --pairwise_del --seed 3 --probs --n_ind 850 --n_sites 231488  --labels $LABEL --out $BASE_DIR/No_header_cunner_merged_subsampled_1out50.dist
#Performs MDS using get_PCA.R:
conda activate r_lcpipe
tail -n +3 $BASE_DIR/No_header_cunner_merged_subsampled_1out50.dist | Rscript --vanilla --slave /projects/mjolnir1/people/sjr729/get_PCA.R --no_header --data_symm -n 10 -m mds -o $BASE_DIR/No_header_cunner_merged_subsampled_1out50.mds
```


# 4. Heterozygosity with BEAGLE as input
```bash
module load angsd 

BASE_DIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Cunner"
BEAGLE="$BASE_DIR/cunner_merged.beagle.gz"
GENOME="/projects/mjolnir1/people/sjr729/Skmer_ms/Cunner/genome_assemblies/GCA_020745685.1_fTauAds1.pri.cur_genomic.fna"


## Format the genome
samtools faidx $GENOME
picard CreateSequenceDictionary \
    R=/projects/mjolnir1/people/sjr729/Skmer_ms/Cunner/genome_assemblies/GCA_020745685.1_fTauAds1.pri.cur_genomic.fna \
    O=/projects/mjolnir1/people/sjr729/Skmer_ms/Cunner/genome_assemblies/GCA_020745685.1_fTauAds1.pri.cur_genomic.dict
bwa index $GENOME
## Get saf file
# Define the number of individuals
N_IND=850
# Run angsd without the -GL option
angsd -beagle $BEAGLE \
      -out $BASE_DIR/Cunner_oct24 \
      -doSaf 1 \
      -anc $GENOME \
      -ref $GENOME \
      -nThreads 10 \
      -fai ${GENOME}.fai
      
## Get SFS from saf

## Generate per site thetas


#For Oedulis 
rsync -av --files-from=/home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/bamlist_EUostrea.txt username@remote.cluster:/ /local/destination/path/
