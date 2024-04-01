### Establish Angsd filters values for -minInd -setMinDepthInd -setMinDepth for Clupea


## Per position depth summed across all individuals

```bash
module load angsd
REF="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
BAMLIST="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/01_infofiles/Atmore_modern_bam_list.txt"  # Updated bamlist variable
BASEDIR=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/
OUTPUTFOLDER="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/SetFilters"

angsd \
-bam $BAMLIST \
-ref $REF \
-out $OUTPUTFOLDER/Mar24_minq20_minmapq20 \
-GL 1 -doMajorMinor 1 -doMaf 1 \
-doDepth 1 -doCounts 1 -maxDepth 6000 -dumpCounts 1 \
-minMapQ 20 -minQ 20 \
-remove_bads 1 -only_proper_pairs 1 \
-nThreads 40
```

Plot the depth distribution obtained from ANGSD to establish depth filter with Rscripts/SetAngsdFilters.R

## Count depth 