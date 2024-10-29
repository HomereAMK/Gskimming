
# 1. ngsdist to get p-distance from ANGSD beagle after SNP "calling".
# 1.1 Oyster
```bash


#Get the label list from the bam list
awk '{split($0,a,"/"); split(a[12],b,"_"); print b[1]}' /projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2/docs/bamlist.txt > /projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2/docs/bamlist.labels

BEAGLE=/projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2/angsd/snp_calling_global/combined.subsetted.beagle.gz
LABEL=/projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2/docs/bamlist.labels

#ngsdit
ngsDist --n_threads 40 --geno $BEAGLE --pairwise_del --seed 3 --probs --n_ind 93 --n_sites 3159635 --labels $LABEL --out /projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2/angsd/snp_calling_global/combined.subsetted.dist


#Performs MDS using get_PCA.R:
conda activate r_lcpipe
tail -n +3 /projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2/angsd/snp_calling_global/combined.subsetted.dist | Rscript --vanilla --slave /projects/mjolnir1/people/sjr729/get_PCA.R --no_header --data_symm -n 10 -m mds -o /projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2/angsd/snp_calling_global/combined.subsetted.mds
```
# 1.2 Herring

# 2. filtered out of probable inversion on chr6 beagle file
```bash
# Decompress the original Beagle file
zcat angsd/snp_calling_global/combined.subsetted.beagle.gz > combined.subsetted.beagle

# Filter out the SNPs from chromosome 6
grep -v "NC_079169.1" combined.subsetted.beagle > combined.subsetted.filtered.beagle
# Filter out the SNPs from chromosome 4, 5, 6, 8
grep -v -e "NC_079167.1" -e "NC_079168.1" -e "NC_079169.1" -e "NC_079171.1" combined.subsetted.beagle > combined.subsetted.filtered_chr4568.beagle

# Recompress the filtered file
gzip combined.subsetted.filtered.beagle
gzip angsd/snp_calling_global/combined.subsetted.filtered_chr4568.beagle

#variables
BEAGLE=/projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2/angsd/snp_calling_global/combined.subsetted.filtered_chr4568.beagle.gz
zcat $BEAGLE | tail -n +2 | wc -l

LABEL=/projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2/docs/bamlist.labels

#2825020
#1872976
#ngsdit
ngsDist --n_threads 40 --geno $BEAGLE --pairwise_del --seed 3 --probs --n_ind 93 --n_sites 1872976 --labels $LABEL --out /projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2/angsd/snp_calling_global/combined.subsetted.filtered_chr4568.dist


#Performs MDS using get_PCA.R:
conda activate r_lcpipe
tail -n +3 /projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2/angsd/snp_calling_global/combined.subsetted.filtered_chr4568.dist | Rscript --vanilla --slave /projects/mjolnir1/people/sjr729/get_PCA.R --no_header --data_symm -n 10 -m mds -o /projects/mjolnir1/people/sjr729/SteinOedulis_loco_run2/angsd/snp_calling_global/combined.subsetted.filtered_chr4568.mds