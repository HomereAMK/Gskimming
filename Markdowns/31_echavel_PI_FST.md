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
for query in ${POP[*]}
    do

        angsd -nThreads 12 -ref $REF -anc $REF -bam $DIR/bamlists/${query}_bee_popgen.list -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -setMinDepthInd 1 -setMinDepth 146 -setMaxDepth 400  -isHap 0 -GL 1 -doSaf 1 -out $DATAOUTPUT/${query}_146--400_isHap
    done
done


#by COLONY
OUTPUTDIR=/projects/tomg/people/sjr729/echarvel_bam/bee/4x_PopGen/doSaf_colony
mkdir -p $OUTPUTDIR
for query in ${COLONY[*]}
    do

        angsd -nThreads 12 -ref $REF -anc $REF -bam $DIR/bamlists/${query}_bee_popgen.list -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -setMinDepthInd 1 -setMinDepth 146 -setMaxDepth 400  -isHap 0 -GL 1 -doSaf 1 -out $DATAOUTPUT/${query}_146--400_isHap
    done
done

```