# Gskimming
 Analyses of Nuclear Reads Obtained using Genome skimming

# Module list in Mjolnir 
```  1) skmer/3.3.0   2) apples/20231016   3) fasttree/2.1.11   4) mash/2.3   5) fastme/2.1.6.1  /projects/mjolnir1/people/sjr729/softwares/guppy-64 ```

#if CONSULT used, it needs a lot of memory
```
#salloc --qos=debug --nodes=1 -c 20 --mem-per-cpu 18000 -t 140000
srun --cpus-per-task=10 --mem=18G --time=10:00:00 --pty bash


module load sratoolkit/3.0.0  sra-tools/3.0.3
```

# BBmap from Skmer pipeline try for the Phacochoerus dataset (fastq not gzip)
```
#Make the fastq list
ls /projects/mjolnir1/people/sjr729/Phacochoerus/*.fastq | cut -d '_' -f 1 > Phacochoerus_list.txt
```

```
#!/bin/bash
#SBATCH --job-name=BBmap_Phacochoerus
#SBATCH --output=BBmap-Phacochoerus.out
#SBATCH --error=BBmap-Phacochoerus.err
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=180G    # memory per cpu-core
#SBATCH --time=05:00:00
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=homerejalves.monteiro@sund.ku.dk

#activate proper conda env
conda activate tuto

for x in `cat /projects/mjolnir1/people/sjr729/Phacochoerus/Phacochoerus_list.txt`; do 
        /projects/mjolnir1/people/sjr729/softwares/skimming_scripts/bbmap_pipeline.sh ${x}_1.fastq ${x}_2.fastq ${x}_merged.fastq
done 

```

accession_list=(
...
ERR9709216
ERR9709220
ERR9709222
...
)

```

for accession in "${accession_list[@]}"; do
    prefetch "$accession"
    fasterq-dump --split-files --outdir /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Clupea/ "$accession"
done
```
#gzip all the fastq
find /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Phacochoerus -type f -name "*.fastq" | xargs -I {} -P 30 gzip -k {}


```
#test on few clupea
#(base) [sjr729@mjolnirhead01fl skimming_scripts]$ pwd
/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts


#Linking some fastq to a new test dir 
mkdir testPhacochoerus

test_Phacochoeruslist=(
SRR19174500_1.fastq.gz
SRR19174500_2.fastq.gz
SRR19174498_1.fastq.gz
SRR19174498_2.fastq.gz
SRR19174497_1.fastq.gz
SRR19174497_2.fastq.gz
SRR19174502_1.fastq.gz
SRR19174502_2.fastq.gz
)
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testPhacochoerus
for file in "${test_Phacochoeruslist[@]}"; do
    ln -s /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Phacochoerus/"$file" .
done

#Performing 
conda activate tutorial
bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testPhacochoerus
```



```
#Clupea skimming processing
srun --cpus-per-task=30 --mem=180G --time=96:00:00 --pty bash

cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
ln -s /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Clupea/*.fastq.gz ./
conda activate tutorial
bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testPhacochoerus -r 30 -f 30
```

#skims_processing_pipeline.sh
#Mandatory inputs:
#-x  path to folder containing reads (split reads that need to be merged)
#-g  path to reference library
#-r  threads for RESPECT; default: 8
#-d  number of iteration cycles for RESPECT, default: 1000
#-f  number of cores for SKMER, default: 8"