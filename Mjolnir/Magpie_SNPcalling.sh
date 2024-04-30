

## Magpie set angsd filters


```bash
REF="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/angsd/01_infofiles/GCA_013398635.1_ASM1339863v1_genomic.fna"
BAMLIST="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/angsd/01_infofiles/Magpie_bam_list.txt"  # Updated bamlist variable
BASEDIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/angsd/"
OUTPUTFOLDER="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/angsd/SetFilters"  

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

## Magpie SNP calling
```bash
SEARCH_DIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/angsd/realigned/"
LIST_FILE="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/Atmore_modern_full_paths.txt"
OUTPUT_FILE="Atmore_modern_bam_list.txt"
BAMS=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/angsd/realigned/*.bam
ls $BAMS > 01_infofiles/Magpie_bam_list.txt

GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/angsd/01_infofiles/GCA_013398635.1_ASM1339863v1_genomic.fna"
BAMLIST="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/angsd/01_infofiles/Magpie_bam_list.txt"  # Updated bamlist variable
BASEDIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/angsd/"
OUTPUTFOLDER="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/angsd/SNP_calling"  # Updated output folder variable
N_IND=$(cat $BAMLIST | wc -l)  # Count the number of individuals based on the bamlist

# SNP Variant Calling with ANGSD
angsd \
-bam $BAMLIST \
-ref $GENOME \
-out "${OUTPUTFOLDER}/mar24_Magpie_SNP_minInd0.25_MinDepth450_MaxDepth1600" \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
-minMapQ 20 -minQ 20 -minInd $((N_IND*1/4)) \
-setMinDepthInd 1 -setMinDepth 450 -setMaxDepth 1600 \
-doCounts 1 -dumpCounts 2 \
-GL 1 -doGlf 2 \
-doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 -rmTriallelic 0.05 -doPost 1 -doGeno 8 \
-doIBS 1 -doCov 1 -makeMatrix 1 \
-nThreads 40
```

