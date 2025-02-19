



## Get lists of samples
```bash
#!/bin/bash

# Array of individuals
IND=("CLEW_03" "HYPP_04" "HYPP_06" "NISS_08" "NISS_13" "ORIS_17" "ORIS_19" "PONT_03" "PONT_13" "TRAL_19" "VAGS_06" "VAGS_13" "WADD_12" "WADD_13")

# Paths
BAMLIST=/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/locopipe/docs/bamlist.txt
REF=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/Tanguy/fileOegenome10scaffoldC3G.fasta
OUTDIR="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/HC/PopGenEstimates"

# Read BAM filenames into an array
mapfile -t BAMS < "$BAMLIST"

for i in "${!IND[@]}"; do
    ind=${IND[$i]}
    bam=${BAMS[$i]}

    echo "Running ANGSD for $ind with BAM file $bam..."

    angsd -nThreads 30 \
          -ref "$REF" \
          -anc "$REF" \
          -i "$bam" \
          -uniqueOnly 1 \
          -remove_bads 1 \
          -only_proper_pairs 1 \
          -trim 0 \
          -C 50 \
          -baq 1 \
          -minMapQ 30 \
          -minQ 30 \
          -setMinDepthInd 1 \
          -GL 1 \
          -doSaf 1 \
          -out "$OUTDIR/$ind"
done

```
## Runs realSFS
```bash
# Loop through each individual again and compute SFS from the generated SAF files
for i in "${!IND[@]}"; do
    ind=${IND[$i]}
    
    echo "Computing SFS for $ind using $OUTDIR/$ind.saf.idx..."
    
    # Run realSFS on each individual's SAF index
    realSFS -P 40 -fold 1 "$OUTDIR/$ind.saf.idx" > "$OUTDIR/$ind.sfs"
    
done

echo "All SFS computations complete."
```

## Calculates the thetas for each site:
```bash
# Loop through each individual and generate PopGenEstimates
for i in "${!IND[@]}"; do
    ind=${IND[$i]}
    
    echo "Running saf2theta for $ind using $OUTDIR/$ind.saf.idx and $OUTDIR/$ind.sfs..."
    
    realSFS saf2theta \
        "$OUTDIR/$ind.saf.idx" \
        -sfs "$OUTDIR/$ind.sfs" \
        -fold 1 \
        -outname "$OUTDIR/$ind.PopGenEstimates"
    
done

echo "All PopGenEstimates computations complete."
```


## Gets summary of thetas:
```bash
# Loop through each individual and print the thetas summary
#for i in "${!IND[@]}"; do
#    ind=${IND[$i]}
    
#    echo "Printing thetas summary for $ind..."
#    thetaStat print "$OUTDIR/$ind.PopGenEstimates.thetas.idx" > "$OUTDIR/$ind.PopGenEstimates.Print"
#done

#echo "All thetas summaries printed."
```

## Performs final calculations:
```bash
#for i in "${!IND[@]}"; do
#    ind=${IND[$i]}
#    bam=${BAMS[$i]}

#    echo "Processing $ind..."
#Call R script, passing in the .Print file, N_IND, and the ind name:
    #    (update the path/filenames to match your actual usage)
#    Rscript --vanilla --slave \
#      /maps/projects/mjolnir1/people/sjr729/Skmer_ms/scripts/GetsThetaSummaries.R \
#      "$OUTDIR/${ind}.PopGenEstimates.Print" \
#      1 \
#      "$ind"

# Redirect **all** stdout from the loop into a single file:
#done > "$OUTDIR/Dec24--HC_Oedulis.PopGenEstimates.txt"
```

## Loop using thetaStat do_stat 
```bash
for i in "${!IND[@]}"; do
    ind=${IND[$i]}

    echo "Calculating average theta for $ind..."

    thetaStat do_stat \
      "$OUTDIR/$ind.PopGenEstimates.thetas.idx" \
      -nChr 2 \
      > "$OUTDIR/$ind.PopGenEstimates.AvgTheta"

    echo "Calculating window-based theta for $ind..."

    thetaStat do_stat \
      "$OUTDIR/$ind.PopGenEstimates.thetas.idx" \
      -nChr 2 \
      -win 50000 \
      -step 10000 \
      -outnames "$OUTDIR/$ind.PopGenEstimates.Window"

    # end up with these files:
    # $OUTDIR/$ind.PopGenEstimates.AvgTheta
    # $OUTDIR/$ind.PopGenEstimates.Window.pestPG

done
```











## Gets thetas per sliding window:
```bash
for query in ${POP[*]}

do
    thetaStat do_stat /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Ind581_${query}_PopGenEstimates.thetas.idx -win 20000 -step 20000 -outnames /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Ind581__20K_${query}_PopGenEstimates-Windows
done
```

## Edits thetas per sliding window:
```bash
for query in ${POP[*]}

do
    cut -f 2,3,4,5,9,14 /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Ind581__20K_${query}_PopGenEstimates-Windows.pestPG | tail -n +2 | sed -r 's/LG//g' | sed 's/^0*//' | awk '$6 > 0' | awk '{print $1"\t"$1":"$2"\t"$2-20000"\t"$2"\t"$6"\t"$3"\t"$3/$6"\t"$4"\t"$4/$6"\t"$5}' | awk 'BEGIN{print "CHR\tSNP\tgPoint\tWindow\tNumberOfSites\tsumTw\tTw\tsumTp\tTp\tTd"}1' > /home/projects/dp_00007/data/hmon/angsd_PopGen/Dec22--Ind581__20K_${query}_PopGenEstimates-Windows.tsv
done
```




/maps/projects/mjolnir1/people/sjr729/Skmer_ms/scripts/