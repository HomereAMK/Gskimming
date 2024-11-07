# 1. Cunner
```bash
SIAVASH_NGSDIST=/projects/mjolnir1/people/sjr729/Skmer_ms/scripts/ngsDist/ngsDist
BASE_DIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Cunner/ngsdist"
BEAGLE="$BASE_DIR/cunner_merged.beagle.gz"
LABEL="/projects/mjolnir1/people/sjr729/Skmer_ms/Cunner/ngsdist/cunner_labels.tsv"

#downsampling
#Prepare a geno file by subsampling one SNP in every 20 SNPs in the beagle file
#zcat $BEAGLE | awk 'NR % 20 == 0' | cut -f 4- | gzip  > $BASE_DIR/cunner_merged_subsampled_1out20.beagle.gz
zcat $BEAGLE | awk 'NR % 20 == 0' | gzip > $BASE_DIR/cunner_merged_subsampled_1out20.beagle.gz


#number of sites
DOWNSAMP_BEAGLE="$BASE_DIR/cunner_merged_subsampled_1out20.beagle.gz"
zcat $DOWNSAMP_BEAGLE | tail -n +2 | wc -l
#DOWNSAMP_BEAGLE="$BASE_DIR/cunner_merged_subsampled_1out50.beagle.gz"

#ngsdist
$SIAVASH_NGSDIST --n_threads 30 --geno $DOWNSAMP_BEAGLE --pairwise_del --avg_pi --evol_model 0 --seed 3 --probs --n_ind 850 --n_sites 578721 --labels $LABEL --out $BASE_DIR/Siavash_ngsdist_cunner_merged_subsampled_1out20_avg_pi_evol_model_0.dist

#ERROR: [read_geno] GENO file not at EOF. Check GENO file and number of sites!
#        : Numerical result out of range

$SIAVASH_NGSDIST --n_threads 30 --geno $DOWNSAMP_BEAGLE --pairwise_del --avg_pi --evol_model 0 --seed 3 --probs --n_ind 850 --n_sites 578721 --labels $LABEL --out $BASE_DIR/Siavash_ngsdist_cunner_merged_subsampled_1out20_avg_pi_evol_model_0.dist

$SIAVASH_NGSDIST \
--n_threads 10 --geno $DOWNSAMP_BEAGLE \
--seed 3 --probs --n_ind 850 \
--n_sites 578721 --labels $LABEL  \
--theta --tot_sites 732626995 \
--out $BASE_DIR/Siavash_ngsdist_cunner_merged_subsampled_1out20_totsites_theta.dist 


zcat $BEAGLE | awk 'NR % 10 == 1' | gzip > $BASE_DIR/cunner_1out10.beagle.gz
DOWNSAMP_BEAGLE="$BASE_DIR/cunner_1out10.beagle.gz"
zcat $DOWNSAMP_BEAGLE | tail -n +2 | wc -l #1157443
#done
#N_SITES=$(( (732626995/10)/(11574435/1157443) ))
whole_genome_sites=732626995
N_SITES=$(( whole_genome_sites / 10 ))


$SIAVASH_NGSDIST \
--n_threads 10 --geno $DOWNSAMP_BEAGLE \
--seed 3 --probs --n_ind 850 \
--n_sites 1157443 --labels $LABEL \
--theta --evol_model 0 --tot_sites $N_SITES  \
--out $BASE_DIR/Siavash_ngsdist_cunner_1out10_theta_--evol_model0_nsites$N_SITES.dist



```

# 2. Oyster
```bash
SIAVASH_NGSDIST=/projects/mjolnir1/people/sjr729/Skmer_ms/scripts/ngsDist/ngsDist
BASE_DIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/ngsdist"
BEAGLE="$BASE_DIR/30jan23_prunedLDminweight0.5_PopStruct.beagle.gz"
LABEL="$BASE_DIR/bamlist_EUostrea.labels"
BEAGLE_NOHEADER="30jan23_prunedLDminweight0.5_PopStruct.beagle_noheader.gz"


#$SIAVASH_NGSDIST --n_threads 30 --geno $BEAGLE --pairwise_del --seed 3 --probs --n_ind 581 --n_sites 1404180 --labels $LABEL --out $BASE_DIR/Siavash_Oedulis_beagle_avg_pi_evol_model_0.dist --avg_pi --evol_model 0

#$SIAVASH_NGSDIST --n_threads 30 --geno $BEAGLE --pairwise_del --avg_pi --evol_model 0 --seed 3 --probs --n_ind 581 --n_sites 1404181 --labels $LABEL --out $BASE_DIR/Siavash_ngsdist_oedulis_merged_subsampled_1out20_avg_pi_evol_model_0.dist



#downsampling
zcat $BEAGLE | awk 'NR % 20 == 0' | gzip > $BASE_DIR/oedulis_1out20.beagle.gz
DOWNSAMP_BEAGLE="$BASE_DIR/oedulis_1out20.beagle.gz"
zcat $DOWNSAMP_BEAGLE | tail -n +2 | wc -l
zcat $BEAGLE | tail -n +2 | wc -l #1140180
zcat $BEAGLE_NOHEADER | tail -n +2 | wc -l

#done
$SIAVASH_NGSDIST --n_threads 30 --geno $DOWNSAMP_BEAGLE --pairwise_del --avg_pi --evol_model 0 --seed 3 --probs --n_ind 581 --n_sites 57009 --labels $LABEL --out $BASE_DIR/Siavash_ngsdist_oedulis_merged_subsampled_1out20_avg_pi_evol_model_0.dist



$SIAVASH_NGSDIST \
--n_threads 10 --geno $DOWNSAMP_BEAGLE \
--seed 3 --probs --n_ind 581 \
--n_sites 57009 --labels $LABEL  \
--theta --tot_sites 840354420 \
--out $BASE_DIR/Siavash_ngsdist_oedulis_merged_subsampled_1out20_theta_totsites.dist

$SIAVASH_NGSDIST \
--n_threads 10 --geno $DOWNSAMP_BEAGLE \
--seed 3 --probs --n_ind 581 \
--n_sites 57009 --labels $LABEL  \
--theta --pairwise_del \
--out $BASE_DIR/Siavash_ngsdist_oedulis_merged_subsampled_1out20_theta_pairwisedel.dist





zcat $BEAGLE | awk 'NR % 10 == 1' | gzip > $BASE_DIR/oedulis_1out10.beagle.gz
DOWNSAMP_BEAGLE="$BASE_DIR/oedulis_1out10.beagle.gz"
zcat $DOWNSAMP_BEAGLE | tail -n +2 | wc -l
#done
#N_SITES=$(( (840354420/10)/(5684643/1140180) ))
N_SITES=$(( 840354420*1140180/10/5684643 ))


$SIAVASH_NGSDIST \
--n_threads 10 --geno $DOWNSAMP_BEAGLE \
--seed 3 --probs --n_ind 581 \
--n_sites 114019 --labels $LABEL \
--theta --evol_model 0 --tot_sites $N_SITES  \
--out $BASE_DIR/Siavash_ngsdist_oedulis_merged_subsampled_1out10_theta_--evol_model0_nsites$N_SITES.dist
```             




# 3. Oyster with 14M variavle sites
```bash 
BASE_DIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/ngsdist/varSites_sensitivity"
BEAGLE="$BASE_DIR/Oedulis_varSites_sensitivity_-SNP_pval1e-03_setMinDepth500_setMaxDepth1500.beagle.gz"
BAMLIST="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/ngsdist/bam_molecol_list.txt"
LABEL="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/ngsdist/bam_molecol_list.labels"
SIAVASH_NGSDIST="/projects/mjolnir1/people/sjr729/Skmer_ms/scripts/ngsDist/ngsDist"
OUT_DIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/ngsdist"

#downsamp
zcat $BEAGLE | awk 'NR % 20 == 1' | gzip > $BASE_DIR/Oedulis_varSites_sensitivity_-SNP_pval1e-03_setMinDepth500_setMaxDepth1500_1out20.beagle.gz
DOWNSAMP_BEAGLE="$BASE_DIR/Oedulis_varSites_sensitivity_-SNP_pval1e-03_setMinDepth500_setMaxDepth1500_1out20.beagle.gz"
zcat $DOWNSAMP_BEAGLE | tail -n +2 | wc -l #698453

#labels
sed 's|/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/ngsdist/bam_molecol/||' /projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/ngsdist/bam_molecol_list.txt | cut -c 1-7 > /projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/ngsdist/bam_molecol_list.labels

N_SITES=$((13969066/20))

$SIAVASH_NGSDIST \
--n_threads 10 --geno $DOWNSAMP_BEAGLE \
--seed 3 --probs --n_ind 519 \
--n_sites 698453 --labels $LABEL \
--theta --evol_model 0 --tot_sites $N_SITES  \
--out $OUT_DIR/Siavash_ngsdist_Oedulis_varSites_sensitivity_-SNP_pval1e-03_setMinDepth500_setMaxDepth1500_1out20_theta_--evol_model0_nsites$N_SITES.dist
