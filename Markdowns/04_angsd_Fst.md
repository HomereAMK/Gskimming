### FST matrix with ANGSD pipeline

## Extracting `chr` and `pos` Columns for -sites (ANGSD): Inversion-free SNP list

```sh
zcat SNP_calling/screenfeb24_Atmore_modern_SNP_Inv_free_regions_minInd0.25.pos.gz | cut -f 1,2 > 01_infofiles/sites_ClupeaAtmore_InvFree_file.txt
```

## Creating bamlists for each pop for pairwise FST calculations

```sh
# Define the path to your input file
BAMLIST="01_infofiles/bam_ModernClupea_annot.tsv"

# Define an array of populations
POP=("BLEK" "IDEI" "IDEO" "KALI" "KALM" "KARM" "MORE" "MASE" "RISO")  # Add more population codes as needed

# Loop through each population
for query in "${POP[@]}"
do 
    # Use grep to find lines matching the population code and cut to select only the bam file paths
    grep -w "${query}" "$BAMLIST" | cut -f1 > "01_infofiles/Jan24--ModernClupea_${query}-Fst.list"
done
```

## angsd index sites files that contain info on Maj an Min allele
```sh

gunzip -c SNP_calling/screenfeb24_Atmore_modern_SNP_Inv_free_regions_minInd0.25.mafs.gz | cut -f 1,2,3,4 | tail -n +2 > 01_infofiles/sites_ClupeaAtmore_InvFree_majmin_file.txt
angsd sites index 01_infofiles/sites_ClupeaAtmore_InvFree_majmin_file.txt
```

## Get minor allele frequency estimation from angsd for each population / group
```sh
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
BASEDIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/"
OUTPUTFOLDER="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/FST_based"  # Updated output folder variable
EXTRA_ARG='-remove_bads 1 -only_proper_pairs 1 -C 50'


for i1 in `seq 0 $((${#POP[@]}-2))`
do
    N_IND=`cat 01_infofiles/Jan24--ModernClupea_${POP[i1]}-Fst.list | wc -l`
    angsd -nThreads 10 \
    -bam 01_infofiles/Jan24--ModernClupea_${POP[i1]}-Fst.list \
    -anc $GENOME \
    -ref $GENOME \
    -out FST_based/Feb24--mindInd0.25_Unfolded_ModerClupea_InvFreeSNPlist_${POP[i1]} \
    -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1 \
    -P $THREADS \
    -minInd $((N_IND*1/4)) -minQ 20 -minMapQ 20 \
    -sites 01_infofiles/sites_ClupeaAtmore_InvFree_majmin_file.txt \
    $EXTRA_ARG
done
```

## Estimate Fst in each pair of populations

```sh

OUTPUTFOLDER="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/FST_based"  # Updated output folder variable

cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/FST_based
for i1 in `seq 0 $((${#POP[@]}-2))`
do
    for i2 in `seq $((i1+1)) $((${#POP[@]}-1))`
    do
        pop1="Feb24--mindInd0.25_Unfolded_ModerClupea_InvFreeSNPlist_${POP[i1]}"
        pop2="Feb24--mindInd0.25_Unfolded_ModerClupea_InvFreeSNPlist_${POP[i2]}"
        N_SITES=`realSFS print $pop1.saf.idx $pop2.saf.idx | wc -l`
         echo -ne "${POP[i1]}\t${POP[i2]}\t$N_SITES\t"
        if [[ $N_SITES == 0 ]]; then
            echo "NA"
        else
            realSFS $pop1.saf.idx $pop2.saf.idx -fold 1 -P 40 > Feb24--mindInd0.25_Unfolded_ModerClupea_InvFreeSNPlist_${POP[i1]}.${POP[i2]}.sfs
            realSFS fst index $pop1.saf.idx $pop2.saf.idx -sfs Feb24--mindInd0.25_Unfolded_ModerClupea_InvFreeSNPlist_${POP[i1]}.${POP[i2]}.sfs -fold 1 -P 40 -fstout Feb24--mindInd0.25_Unfolded_ModerClupea_InvFreeSNPlist_${POP[i1]}.${POP[i2]}
            realSFS fst stats Feb24--mindInd0.25_Unfolded_ModerClupea_InvFreeSNPlist_${POP[i1]}.${POP[i2]}.fst.idx -P 40
            fi
        done
done > Mar24--mindInd0.25_Unfolded_ModerClupea_InvFreeSNPlist--Fst.tsv
```
🫡