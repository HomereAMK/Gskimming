## Skmer Preprocessing Mapping and Realignment for AtmoreClupea


```bash
module load angsd parallel
BAMLIST="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/angsd/01_infofiles/Atmore_modern_bam_list.txt"
POP=("BLEK" "IDEI" "IDEO" "KALI" "KALM" "KARM" "MORE" "MASE" "RISO")  # Add more population 
REF="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
OUTPUTDIR="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/angsd/Het/"
```

### Runs ANGSD under -doSaf on all populations, estimate the site allele frequency likelihood for each pop
```bash
for query in ${POP[*]}
    do

        angsd -nThreads 40 -ref $REF -anc $REF -bam /home/projects/dp_00007/people/hmon/EUostrea/01_infofiles/EUostrea_${query}-Fst.list -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -setMinDepthInd 1 -setMinDepth 600 -setMaxDepth 1200 -GL 1 -doSaf 1 -out $OUTPUTDIR/Mar24--Unfolded_Magpie_NIC_${query}
    done
done