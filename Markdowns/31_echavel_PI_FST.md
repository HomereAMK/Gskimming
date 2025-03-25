# BEE PI per colony

## 1. Create bamlists by population and colony
```bash
DIR=/projects/tomg/people/sjr729/echarvel_bam/bee/4x_PopGen
#by pop
awk -F'\t' 'NR>1 {print $2 > ("bamlists/" $3 "_bee_popgen.list")}' bee_pop_col.tsv
#by colony
awk -F'\t' 'NR>1 {print $2 > ("bamlists/" $4 "_bee_popgen.list")}' bee_pop_col.tsv
```

## 2. ANGSD -doSaf per pop and per colony
```bash
module load angsd parallel


POP=("BM" "FAS" "MAS")
COLONY=("C_3145" "C_3152" "C_3154" "C_3155" "C_3156" "C_3160" "C_3162" "C_3163" "C_3164" "C_3165" "C_3169" "C_3172" "C_3173" "C_3175" "C_3498" "C_3504" "C_3505" "C_3507" "C_3508" "C_3510" "C_3516" "C_3517" "C_3518" "C_3519" "C_3520" "C_3521" "C_3522" "C_3523" "C_3525" "C_3780" "C_3800" "C_3803" "C_3822" "C_3824" "C_3828" "C_3841" "C_3858" "C_3862" "C_3867" "C_3870")

REF="/projects/tomg/people/sjr729/echarvel_bam/bee/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
DIR=/projects/tomg/people/sjr729/echarvel_bam/bee/4x_PopGen

#by POP
OUTPUTDIR=/projects/tomg/people/sjr729/echarvel_bam/bee/4x_PopGen/doSaf_pop
mkdir -p $OUTPUTDIR
for query in "${POP[@]}"; do
    angsd -nThreads 10 \
        -ref "$REF" \
        -anc "$REF" \
        -bam "$DIR/bamlists/${query}_bee_popgen.list" \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -only_proper_pairs 1 \
        -trim 0 \
        -C 50 \
        -baq 1 \
        -minMapQ 20 \
        -minQ 20 \
        -setMinDepthInd 1 \
        -isHap 0 \
        -GL 1 \
        -doSaf 1 \
        -out "$OUTPUTDIR/${query}_isHap"
        #-setMinDepth 146 \
        #-setMaxDepth 400 
done

#by COLONY
OUTPUTDIR="/projects/tomg/people/sjr729/echarvel_bam/bee/4x_PopGen/doSaf_colony"
mkdir -p "$OUTPUTDIR"

# Ensure COLONY array is defined (e.g., COLONY=("BM" "FAS" "MAS"))
for query in "${COLONY[@]}"; do
    angsd -nThreads 11 \
        -ref "$REF" \
        -anc "$REF" \
        -bam "$DIR/bamlists/${query}_bee_popgen.list" \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -only_proper_pairs 1 \
        -trim 0 \
        -C 50 \
        -baq 1 \
        -minMapQ 20 \
        -minQ 20 \
        -setMinDepthInd 1 \
        -isHap 0 \
        -GL 1 \
        -doSaf 1 \
        -out "$OUTPUTDIR/${query}_isHap"
        #        -setMinDepth 146 \
#        -setMaxDepth 400 \
done
```

## 3. Run realSFS
```bash
#by POP
POP=("BM" "FAS" "MAS")
REF="/projects/tomg/people/sjr729/echarvel_bam/bee/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
DIR=/projects/tomg/people/sjr729/echarvel_bam/bee/4x_PopGen
INPUTDIR=/projects/tomg/people/sjr729/echarvel_bam/bee/4x_PopGen/doSaf_pop
for query in "${POP[@]}"; 
    do
        realSFS -P 40 -fold 1 $INPUTDIR/${query}_isHap.saf.idx > $INPUTDIR/${query}_isHap.sfs
done


#by colony
COLONY=("C_3145" "C_3152" "C_3154" "C_3155" "C_3156" "C_3160" "C_3162" "C_3163" "C_3164" "C_3165" "C_3169" "C_3172" "C_3173" "C_3175" "C_3498" "C_3504" "C_3505" "C_3507" "C_3508" "C_3510" "C_3516" "C_3517" "C_3518" "C_3519" "C_3520" "C_3521" "C_3522" "C_3523" "C_3525" "C_3780" "C_3800" "C_3803" "C_3822" "C_3824" "C_3828" "C_3841" "C_3858" "C_3862" "C_3867" "C_3870")
REF="/projects/tomg/people/sjr729/echarvel_bam/bee/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
DIR=/projects/tomg/people/sjr729/echarvel_bam/bee/4x_PopGen
INPUTDIR=/projects/tomg/people/sjr729/echarvel_bam/bee/4x_PopGen/doSaf_colony
for query in "${COLONY[@]}"; 
    do
        realSFS -P 40 -fold 1 $INPUTDIR/${query}_isHap.saf.idx > $INPUTDIR/${query}_isHap.sfs
done

```



## 4. Calculates the thetas for each site
```bash
#by POP
POP=("BM" "FAS" "MAS")
REF="/projects/tomg/people/sjr729/echarvel_bam/bee/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
DIR=/projects/tomg/people/sjr729/echarvel_bam/bee/4x_PopGen
INPUTDIR=/projects/tomg/people/sjr729/echarvel_bam/bee/4x_PopGen/doSaf_pop
OUTPUTDIR=/projects/tomg/people/sjr729/echarvel_bam/bee/4x_PopGen/doSaf_pop/theta
mkdir -p $OUTPUTDIR

for query in "${POP[@]}"; 
    do
        realSFS saf2theta $INPUTDIR/${query}_isHap.saf.idx -sfs $INPUTDIR/${query}_isHap.sfs -fold 1 -outname $OUTPUTDIR/${query}_isHap_PopGenEstimates
done

#by COLONY
...

```
## 5. Gets summary of thetas:

```bash
# Define paths
DIR="/projects/tomg/people/sjr729/echarvel_bam/bee/4x_PopGen"
POP=("BM" "FAS" "MAS")
INPUTDIR="/projects/tomg/people/sjr729/echarvel_bam/bee/4x_PopGen/doSaf_pop/theta"

# Check input files and permissions
for query in "${POP[@]}"; do
    if [ ! -f "$INPUTDIR/${query}_isHap_PopGenEstimates.thetas.idx" ]; then
        echo "ERROR: Missing $INPUTDIR/${query}_isHap_PopGenEstimates.thetas.idx"
        exit 1
    fi
    if [ ! -f "$DIR/bamlists/${query}_bee_popgen.list" ]; then
        echo "ERROR: Missing $DIR/bamlists/${query}_bee_popgen.list"
        exit 1
    fi
done

if [ ! -w "$INPUTDIR" ]; then
    echo "ERROR: Cannot write to $INPUTDIR"
    exit 1
fi

# Generate theta.Print files
for query in "${POP[@]}"; do
    thetaStat print "$INPUTDIR/${query}_isHap_PopGenEstimates.thetas.idx" > "$INPUTDIR/${query}_isHap_PopGenEstimates_theta.Print"
done

# Run R script for each population
for query in "${POP[@]}"; do
    N_IND=$(wc -l < "$DIR/bamlists/${query}_bee_popgen.list")
    Rscript --vanilla --slave /projects/tomg/people/sjr729/scripts/GetsThetaSummaries.R \
        "$INPUTDIR/${query}_isHap_PopGenEstimates_theta.Print" \
        "$N_IND" \
        "$query"
done > "$INPUTDIR/POP_bee_Mar25--isHap.PopGenEstimates.txt"