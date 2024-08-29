# Git pull KRANK
```bash
cd /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/KRANK
git pull
module load gcc
make
```



# test KRANK with treshold fixed on 1 ind Oedulis 
```bash
conda activate skimming_echarvel
module purge
module load parallel
SCRIPT_DIR=/projects/mjolnir1/people/sjr729/skimming_scripts-echarvel
OUTPUT_DIRECTORY=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/skmer/Krank_test_12.08

LIBRARIES=$"${SCRIPT_DIR}/KRANK/lib_reps_adpt-k29_w35_h13_b16_s8 cxvvv bbv bv h h h hhh{SCRIPT_DIR}/KRANK/pangenome-05-2024-lib_rand_free-k29_w34_h13_b16_s8"
ls ${LIBRARIES}

${SCRIPT_DIR}/KRANK/krank query \
--library-dir ${LIBRARIES} \
--query-file /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/Stein_Mortensen_fqs/unclassified-kra_SRR26416938_.fq \
--max-match-distance 5 \
--total-vote-threshold 0.05 0.2 \
--num-threads 10 \
--output-dir "${OUTPUT_DIRECTORY}/krank_output/krank_reports/"
```

# Rerun the skimming_pipeline.sh with KRANK treshold fixed
```bash
conda activate skimming_echarvel
module purge
module load parallel

#herring
bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline.sh \
-i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/fastq \
-o /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/echarvel_fast-skims_14.08_FIXEDthreshold0.05_0.2 -t 20 -p 10 \
-f _1.fastq.gz \
-r _2.fastq.gz 

#oedulis
conda activate skimming_echarvel
module purge
module load parallel
bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline.sh \
-i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/Stein_raw_fq/ \
-o /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/skmer/echarvel_Steinoedulis_fast-skims_results_12.08_FIXEDthreshold0.05_0.2 -t 18 -p 2 \
-f _1.fastq.gz \
-r _2.fastq.gz \
-c 4,2,1,0.5,0.25

#
sbatch --mem=50g --cpus-per-task=20 --time=52:00:00 --wrap="conda activate skimming_echarvel; module purge; module load parallel; bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline.sh -i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/Stein_raw_fq/ -o /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/skmer/echarvel_Steinoedulis_fast-skims_results_12.08_FIXEDthreshold0.05_0.2 -t 18 -p 2 -f _1.fastq.gz -r _2.fastq.gz"


#magpie
conda activate skimming_echarvel
module purge
module load parallel
bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline.sh \
-i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/done \
-o /projects/mjolnir1/people/sjr729/Skmer_ms/skmer1/echarvel_skimpip_magpie_KrankFixed_14.08 -t 20 -p 2 \
-f _1.fastq.gz \
-r _2.fastq.gz \
-c 4,2,1,0.5,0.25
```


# Non random library in /projects/mjolnir1/people/sjr729/Skmer_ms (subsampling not working)
```bash
#get the non-rand bacterial lib
 cd /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/KRANK
wget https://ter-trees.ucsd.edu/data/krank/wol_v1-lib_reps_adpt-k29_w35_h13_b16_s8.tar.gz
tar -zxf ./wol_v1-lib_reps_adpt-k29_w35_h13_b16_s8.tar.gz

#on herring with increased threshold 0.2 0.2
conda activate skimming_echarvel
module purge
module load parallel
bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline_nonrandolib.sh \
-i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/fastq \
-o /projects/mjolnir1/people/sjr729/Skmer_ms/Herring/skmer1/echarvel_skimpip_16.08_nonrandlib_tresh0.2_0.2 -t 20 -p 2 \
-f _1.fastq.gz \
-r _2.fastq.gz 

#finalize with the KRANK step hashed
conda activate skimming_echarvel
module purge
module load parallel
bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline_withoutKRANK.sh \
-i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/fastq \
-o /projects/mjolnir1/people/sjr729/Skmer_ms/Herring/skmer1/echarvel_skimpip_16.08_nonrandlib_tresh0.2_0.2 -t 20 -p 2 \
-f _1.fastq.gz \
-r _2.fastq.gz 



#
sbatch --mem=100G --time=75:00:00 --cpus-per-task=20 --wrap="conda activate skimming_echarvel; module purge; module load parallel; bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline_nonrandolib.sh -i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/fastq -o /projects/mjolnir1/people/sjr729/Skmer_ms/Herring/skmer1/echarvel_skimpip_16.08_nonrandlib_tresh0.2_0.2 -t 20 -p 10 -f _1.fastq.gz -r _2.fastq.gz"




#magpie
conda activate skimming_echarvel
module purge
module load parallel
bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline_nonrandolib.sh \
-i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/done \
-o /projects/mjolnir1/people/sjr729/Skmer_ms/Magpie/skmer1/echarvel_magpie_skimpip_results_19.08_nonrandlib_tresh0.2_0.2 -t 20 -p 10 \
-f _1.fastq.gz \
-r _2.fastq.gz 


#oedulis + subsampling
conda activate skimming_echarvel
module purge
module load parallel
bash /projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/skimming_pipeline_nonrandolib.sh \
-i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/oedulis_RawFastq_14.08.24 \
-o /projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/skmer1/echarvel_oedulis_skimpip_results_17.08_nonrandlib_tresh0.2_0.2 -t 20 -p 10 \
-f _1.fastq.gz \
-r _2.fastq.gz \
-c 4,2,1,0.5,0.25
#subsampling_and_estimates.sh on oedulis
conda activate skimming_echarvel 
DATE=$(date +%d.%m)
module load parallel
cd /projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/skmer1/echarvel_oedulis_skimpip_results_17.08_nonrandlib_tresh0.2_0.2
INPUTDIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/oedulis_RawFastq_14.08.24"
cd $INPUTDIR
bash /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/subsample_and_estimate2.sh -i $INPUTDIR -t 16 -c 1,0.5,0.25 2>&1 > subsample2"$DATE"_condaskmer2test-_oedulis_skimpip_results.log
