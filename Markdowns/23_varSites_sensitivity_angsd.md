# 1. make a bamlist
```bash
module load samtools
BAMDIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/ngsdist/bam_molecol"
ls -d "$BAMDIR"/*.bam > bam_molecol_list.txt
find "$BAMDIR" -type f -name "*.bam" -exec samtools index {} \;
```

# 2. Angsd Snp calling relaxed filters
```bash
module load angsd 

BAMLIST="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/ngsdist/bam_molecol_list.txt"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/Tanguy/fileOegenome10scaffoldC3G.fasta"
OUTPUTFOLDER="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/ngsdist/varSites_sensitivity"

cd $OUTPUTFOLDER
angsd \
-bam $BAMLIST \
-ref $GENOME \
-out $OUTPUTFOLDER/Oedulis_varSites_sensitivity_-SNP_pval1e-03_setMinDepth500_setMaxDepth1500 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
-minMapQ 20 -minQ 20 -setMinDepthInd 1 -setMinDepth 500 -setMaxDepth 1500 \
-doCounts 1 -dumpCounts 2 \
-GL 1 -doGlf 2 \
-doMajorMinor 1 -doMaf 1 -SNP_pval 1e-3 -doPost 1 -doGeno 8 \
-doIBS 1 -doCov 1 -makeMatrix 1 \
-nThreads 10
```
```txt
-> Sun Nov  3 17:23:05 2024
        -> Arguments and parameters for all analysis are located in .arg file
        -> Total number of sites analyzed: 840002688
        -> Number of sites retained after filtering: 13969066
        [ALL done] cpu-time used =  504099.78 sec
        [ALL done] walltime used =  171420.00 sec
```

# 3. Generate a snps stats file

```bash
module load angsd 

BAMLIST="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/ngsdist/bam_molecol_list.txt"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/Tanguy/fileOegenome10scaffoldC3G.fasta"
OUTPUTFOLDER="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/ngsdist/varSites_sensitivity"

cd $OUTPUTFOLDER
angsd \
-bam $BAMLIST \
-ref $GENOME \
-out $OUTPUTFOLDER/Oedulis_varSites_pval_all_sites \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 \
-minMapQ 20 -minQ 20 \
-doCounts 1 \
-GL 1 \
-doMajorMinor 1 -doMaf 1 \
-doHWE 1 \
-doSnpStat 1 \
-nThreads 10
```

angsd -bam $BAMLIST \
-ref $GENOME \
-doMajorMinor 1 -doMaf 1 \
-gl 1 -snp_pval 1e-1 -P 10 -dosnpstat 1 \
-out $OUTPUTFOLDER/Oedulis_varSites_pval_1e-1_scaffold6 -doHWE 1 -r scaffold6


