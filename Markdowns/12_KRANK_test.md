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
