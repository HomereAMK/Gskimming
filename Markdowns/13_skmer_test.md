
# skmer test with increased sketch size
# after skimming_pipeline.sh (KRANK + bbmap)
# Oedulis


```bash
conda activate skimming_echarvel
module purge
module load parallel

READS="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/skmer/echarvel_Steinoedulis_fast-skims_results_04.08_onlymicrobe/bbmap_reads/"

skmer reference $READS -t -p 10 -l /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/skmer/echarvel_Steinoedulis_fast-skims_results_04.08_onlymicrobe/library_skmer1_sketch50M_12.08 -s 50000000
```


# Skmer2 on herring

