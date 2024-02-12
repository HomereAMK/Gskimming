
### SNP calling

```bash
module load angsd

GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
BAMLIST="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/01_infofiles/Atmore_modern_bam_list.txt"  # Updated bamlist variable
BASEDIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/"
OUTPUTFOLDER="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/SNP_calling"  # Updated output folder variable
N_IND=$(cat $BAMLIST | wc -l)  # Count the number of individuals based on the bamlist

# SNP Variant Calling with ANGSD
angsd \
-bam $BAMLIST \
-ref $GENOME \
-out "${OUTPUTFOLDER}/feb24_Atmore_modern_SNP_minInd0.25" \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
-minMapQ 20 -minQ 20 -minInd $((N_IND*1/4)) \
-doCounts 1 -dumpCounts 2 \
-GL 1 -doGlf 2 \
-doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 -rmTriallelic 0.05 -doPost 1 -doGeno 8 \
-doIBS 1 -doCov 1 -makeMatrix 1 \
-nThreads 40
```

```bash
sbatch --job-name=Atmore_SNP_calling --output=sbatchAtmore_SNP_9feb.out --error=sbatchAtmore_SNP_9feb.err --ntasks=1 --cpus-per-task=40 --mem=180G --time=130:00:00 --mail-type=begin --mail-type=end --mail-type=fail --mail-user=homerejalves.monteiro@sund.ku.dk --wrap="module load angsd && GENOME='/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna' && BAMLIST='/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/01_infofiles/Atmore_modern_bam_list.txt' && OUTPUTFOLDER='/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/SNP_calling' && N_IND=\$(cat \$BAMLIST | wc -l) && angsd -bam \$BAMLIST -ref \$GENOME -out \"\${OUTPUTFOLDER}/batchfeb24_Atmore_modern_SNP_minInd0.25\" -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd \$((N_IND*1/4)) -doCounts 1 -dumpCounts 2 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 -rmTriallelic 0.05 -doPost 1 -doGeno 8 -doIBS 1 -doCov 1 -makeMatrix 1 -nThreads 40"
```

### pcangsd
```bash
SNP_BEAGLE=
pcangsd -b $SNP_BEAGLE \
--selection \
--minMaf 0.05 \
--sites_save \
--threads 40 \
-e 10 \
-o $OUTPUTFOLDER/feb24_Atmore_modern_SNP_selection_minmaf0.05_pcangsd_e10
```

### Without Chromosome with known inversions in herring 
```bash
sbatch --job-name=Atmore_SNP_calling_excl_chr --output=sbatchAtmore_SNP_excl_chr_9feb.out --error=sbatchAtmore_SNP_excl_chr_9feb.err --ntasks=1 --cpus-per-task=40 --mem=180G --time=130:00:00 --mail-type=begin --mail-type=end --mail-type=fail --mail-user=homerejalves.monteiro@sund.ku.dk --wrap="module load angsd && GENOME='/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna' && BAMLIST='/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/01_infofiles/Atmore_modern_bam_list.txt' && OUTPUTFOLDER='/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/SNP_calling_excl_chr' && N_IND=\$(cat \$BAMLIST | wc -l) && angsd -bam \$BAMLIST -ref \$GENOME -out \"\${OUTPUTFOLDER}/screenfeb24_Atmore_modern_SNP_Inv_free_regions_minInd0.25\" -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd \$((N_IND*1/4)) -doCounts 1 -dumpCounts 2 -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 -rmTriallelic 0.05 -doPost 1 -doGeno 8 -doIBS 1 -doCov 1 -makeMatrix 1 -nThreads 40 -rf 01_infofiles/Inv_free_regions_Clupea.txt"
```
