### Popgen estimates for AtmoreClupea


```bash
module load angsd parallel
BAMLIST="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/01_infofiles/Atmore_modern_bam_list.txt"
BASEDIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/"
POP=("BLEK" "IDEI" "IDEO" "KALI" "KALM" "KARM" "MORE" "MASE" "RISO")  # Add more population 
REF="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
OUTPUTDIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/Het/"
```

## Runs ANGSD under -doSaf on all populations, estimate the site allele frequency likelihood for each pop
```bash
for query in ${POP[*]}
    do

        angsd -nThreads 40 -ref $REF -anc $REF -bam $BASEDIR/01_infofiles/Jan24--ModernClupea_${query}-Fst.list -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -setMinDepth 400 -setMaxDepth 800 -GL 1 -doSaf 1 -doCounts 1 -out $OUTPUTDIR/Mar24--Unfolded_Clupea_NIC_MinDepth400-MaxDepth_800${query}
    done
done

```
🤝

## Runs realSFS
```bash

for query in ${POP[*]}

do
    realSFS -P 40 $OUTPUTDIR/Mar24--Unfolded_Clupea_NIC_MinDepth400-MaxDepth_800${query}.saf.idx -tole 1e-08 -maxIter 500 > $OUTPUTDIR/Mar24--Unfolded_Clupea_NIC_MinDepth400-MaxDepth_800${query}.sfs
    
done
```
🤝
## Calculates the thetas for each site and Gets summary of thetas
```bash
for query in ${POP[*]}
do
    realSFS saf2theta $OUTPUTDIR/Mar24--Unfolded_Clupea_NIC_MinDepth400-MaxDepth_800${query}.saf.idx -sfs $OUTPUTDIR/Mar24--Unfolded_Clupea_NIC_MinDepth400-MaxDepth_800${query}.sfs -anc $REF -outname $OUTPUTDIR/Mar24--Unfolded_Clupea_NIC_MinDepth400-MaxDepth_800${query}_PopGenEstimates
    thetaStat print $OUTPUTDIR/Mar24--Unfolded_Clupea_NIC_MinDepth400-MaxDepth_800${query}_PopGenEstimates.thetas.idx > $OUTPUTDIR/Mar24--Unfolded_Clupea_NIC_MinDepth400-MaxDepth_800${query}_PopGenEstimates.Print
done
```
🤝
## Estimate theta
for query in ${POP[*]}
do

    realSFS saf2theta \
    $OUTPUTDIR/Mar24--Unfolded_Clupea_NIC_MinDepth400-MaxDepth_800${query}.saf.idx \
    -outname $OUTPUTDIR/Mar24--Unfolded_Clupea_NIC_MinDepth400-MaxDepth_800${query}_PopGenEstimates \
    -sfs $OUTPUTDIR/Mar24--Unfolded_Clupea_NIC_MinDepth400-MaxDepth_800${query}.sfs \
    -anc $REF \
    -P 40
done

## Gets thetas per sliding window and Edits thetas per sliding window:
```bash

for query in ${POP[*]}
do
    thetaStat do_stat $OUTPUTDIR/Mar24--Unfolded_Clupea_NIC_MinDepth400-MaxDepth_800${query}_PopGenEstimates.thetas.idx -win 20000 -step 20000 -outnames $OUTPUTDIR/Mar24--Unfolded_Clupea_NIC_MinDepth400-MaxDepth_800${query}_PopGenEstimates-20kWindows
    cut -f 2,3,4,5,9,14 $OUTPUTDIR/Mar24--Unfolded_Clupea_NIC_MinDepth400-MaxDepth_800${query}_PopGenEstimates-20kWindows.pestPG | tail -n +2 | sed -r 's/LG//g' | sed 's/^0*//' | awk '$6 > 0' | awk '{print $1"\t"$1":"$2"\t"$2-20000"\t"$2"\t"$6"\t"$3"\t"$3/$6"\t"$4"\t"$4/$6"\t"$5}' | awk 'BEGIN{print "CHR\tSNP\tgPoint\tWindow\tNumberOfSites\tsumTw\tTw\tsumTp\tTp\tTd"}1' > $OUTPUTDIR/Mar24--Unfolded_Clupea_NIC_MinDepth400-MaxDepth_800${query}_PopGenEstimates-20kWindows.tsv
done
```


## heterozygosity 
# get the mafs file with depth filters and Generates a .bed file based on the .mafs file:
```bash
N_IND=`cat /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/01_infofiles/Atmore_modern_bam_list.txt | wc -l`
REF="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
BAMLIST=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/01_infofiles/Atmore_modern_bam_list.txt
OUTPUTDIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/Het/"

angsd -nThreads 40 -ref $REF -bam $BAMLIST -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
-doCounts 1 -dumpCounts 2 \
-minMapQ 20 -minQ 20 -setMinDepth 400 -setMaxDepth 800 \
-GL 1 -doMajorMinor 1 -doMaf 1 \
-out $OUTPUTDIR/Apr24_Clupea_MinDepth400-MaxDepth_800_het
🤝

module load bedtools
zcat $OUTPUTDIR/Apr24_Clupea_MinDepth400-MaxDepth_800_het.mafs.gz | cut -f1,2 | tail -n +2 | awk '{print $1"\t"$2-1"\t"$2}' | bedtools merge -i - > $OUTPUTDIR/Apr24_Clupea_MinDepth400-MaxDepth_800_het.bed
🤝
```

# Creates a position file based on this new .bed:
# Indexs the .pos file created above:
# Runs ANGSD under -doSaf::

```bash
awk '{print $1"\t"($2+1)"\t"$3}'  $OUTPUTDIR/Apr24_Clupea_MinDepth400-MaxDepth_800_het.bed > $OUTPUTDIR/Apr24_Clupea_MinDepth400-MaxDepth_800_het_bed.pos
angsd sites index $OUTPUTDIR/Apr24_Clupea_MinDepth400-MaxDepth_800_het_bed.pos
🤝
wait 
parallel --plus angsd -i {} -ref $REF -anc $REF -sites $OUTPUTDIR/Apr24_Clupea_MinDepth400-MaxDepth_800_het_bed.pos -GL 1 -doSaf 1 -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 20 -minQ 20 -out $OUTPUTDIR/perInd_Het/{/...} :::: $BAMLIST
🤝
```

# Calculates fractions:
# Calculates heterozygosity:
```bash
parallel --tmpdir $OUTPUTDIR/perInd_Het/ -j 1 --plus "realSFS -fold 1 -P 10 {} > $OUTPUTDIR/perInd_Het/{/..}.het" ::: $OUTPUTDIR/perInd_Het/*.saf.idx
🤝
#alternative 
# Run 40 instances of realSFS in parallel using specified tmpdir and plus options 
parallel --tmpdir $OUTPUTDIR/perInd_Het/ --plus "realSFS -fold 1 -P 10 {} > $OUTPUTDIR/perInd_Het/{/..}.het" -j 1 --memfree 160000 ::: $OUTPUTDIR/perInd_Het/*.saf.idx
#change the input!!!!!!!!!! *.het
fgrep '.' *.het | tr ":" " " | awk '{print $1"\t"$3/($2+$3)*100}' | sed -r 's/.het//g' | awk '{split($0,a,"_"); print $1"\t"a[1]"\t"$2"\t"$3'} | awk '{split($2,a,"-"); print $1"\t"a[2]"\t"$3'} > $OUTPUTDIR/Apr24_Clupea_MinDepth400-MaxDepth_8000--HET.Apr24.txt
```




## Heterozygosity Therkildsen lab : noTrans + no chr with known inversions.
```bash
#!/bin/bash
module load angsd parallel

# Variables from your input
BAMLIST="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/01_infofiles/Atmore_modern_bam_list.txt"
BASEDIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/"
POP=("BLEK" "IDEI" "IDEO" "KALI" "KALM" "KARM" "MORE" "MASE" "RISO")  # Populations - customize as needed
REF="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
OUTPUTDIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/Het_therk/"
THREADS=20
MINMAPQ=20
MINQ=33
MINDP=2
MAXDP=800


for LINE in `cat $BAMLIST`; do
    NAME_TEMP=`basename ${LINE%.*}`
    OUTBASE=${NAME_TEMP}'_mindp'$MINDP'_maxdp'$MAXDP'_minq'$MINQ'_minmapq'$MINMAPQ

    ## Get saf file
    angsd \
    -i $LINE \
    -anc $REF \
    -ref $REF \
    -out $OUTPUTDIR$OUTBASE \
    -doSaf 1 \
    -GL 1 \
    -P $THREADS \
    -doCounts 1 -dumpCounts 2 \
    -setMinDepth $MINDP \
    -setMaxDepth $MAXDP \
    -minQ $MINQ \
    -minMapQ $MINMAPQ \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -noTrans 1 -rf /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/01_infofiles/Inv_free_regions_Clupea.txt


    ## Get SFS from saf
    realSFS \
    $OUTPUTDIR$OUTBASE'.saf.idx' \
    -P $THREADS \
    > $OUTPUTDIR$OUTBASE'.ml'

    ## Generate per site thetas
    realSFS saf2theta \
    $OUTPUTDIR$OUTBASE'.saf.idx' \
    -outname $OUTPUTDIR$OUTBASE \
    -sfs $OUTPUTDIR$OUTBASE'.ml' \
    -anc $REF \
    -P $THREADS 

    ## Print per site thetas
    thetaStat print \
    $OUTPUTDIR$OUTBASE'.thetas.idx' | \
    gzip \
    > $OUTPUTDIR$OUTBASE'.thetas.tsv.gz'

    ## Print average thetas
    thetaStat do_stat \
    $OUTPUTDIR$OUTBASE'.thetas.idx' \
    -outnames $OUTPUTDIR$OUTBASE'.average_thetas'
done


```

## Heterozygosity Therkildsen lab : with Trans + no chr with known inversions.
```bash
#!/bin/bash
module load angsd parallel

# Variables from your input
BAMLIST="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/01_infofiles/Atmore_modern_bam_list.txt"
BASEDIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/"
POP=("BLEK" "IDEI" "IDEO" "KALI" "KALM" "KARM" "MORE" "MASE" "RISO")  # Populations - customize as needed
REF="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
OUTPUTDIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/Het_therk_withTrans/"
THREADS=20
MINMAPQ=20
MINQ=33
MINDP=2
MAXDP=800


for LINE in `cat $BAMLIST`; do
    NAME_TEMP=`basename ${LINE%.*}`
    OUTBASE=${NAME_TEMP}'_mindp'$MINDP'_maxdp'$MAXDP'_minq'$MINQ'_minmapq'$MINMAPQ'_withTrans'

    ## Get saf file
    angsd \
    -i $LINE \
    -anc $REF \
    -ref $REF \
    -out $OUTPUTDIR$OUTBASE \
    -doSaf 1 \
    -GL 1 \
    -P $THREADS \
    -doCounts 1 -dumpCounts 2 \
    -setMinDepth $MINDP \
    -setMaxDepth $MAXDP \
    -minQ $MINQ \
    -minMapQ $MINMAPQ \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -rf /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/01_infofiles/Inv_free_regions_Clupea.txt


    ## Get SFS from saf
    realSFS \
    $OUTPUTDIR$OUTBASE'.saf.idx' \
    -P $THREADS \
    > $OUTPUTDIR$OUTBASE'.ml'

    ## Generate per site thetas
    realSFS saf2theta \
    $OUTPUTDIR$OUTBASE'.saf.idx' \
    -outname $OUTPUTDIR$OUTBASE \
    -sfs $OUTPUTDIR$OUTBASE'.ml' \
    -anc $REF \
    -P $THREADS 

    ## Print per site thetas
    thetaStat print \
    $OUTPUTDIR$OUTBASE'.thetas.idx' | \
    gzip \
    > $OUTPUTDIR$OUTBASE'.thetas.tsv.gz'

    ## Print average thetas
    thetaStat do_stat \
    $OUTPUTDIR$OUTBASE'.thetas.idx' \
    -outnames $OUTPUTDIR$OUTBASE'.average_thetas'
done


```

