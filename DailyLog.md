# Gskimming
 Analyses of Nuclear Reads Obtained using Genome skimming

# Module list in Mjolnir 
```  1) skmer/3.3.0   2) apples/20231016   3) fasttree/2.1.11   4) mash/2.3   5) fastme/2.1.6.1  /projects/mjolnir1/people/sjr729/softwares/guppy-64 ```

#if CONSULT used, it needs a lot of memory
```
#salloc --qos=debug --nodes=1 -c 20 --mem-per-cpu 18000 -t 140000
srun --cpus-per-task=20 --mem=80g --time=40:00:00 --pty bash
srun --cpus-per-task=11 --mem=100g --time=120:00:00 --pty bash

srun --cpus-per-task=30 --mem=50g --time=10:00:00 --pty bash


sacct -u sjr729 --units=G --format "JobID%25,JobName%30,Partition,NodeList,Elapsed,CPUTime,ReqMem,MaxRSS,State%15, AllocTRES%32_" --starttime 2024-03-20
module load sratoolkit/3.0.0  sra-tools/3.0.3
```

# BBmap from Skmer pipeline try for the Phacochoerus dataset (fastq not gzip)
```
#Make the fastq list
ls /projects/mjolnir1/people/sjr729/Phacochoerus/*.fastq | cut -d '_' -f 1 > Phacochoerus_list.txt
```

```bash
#!/bin/bash
#SBATCH --job-name=BBmap_Phacochoerus
#SBATCH --output=BBmap-Phacochoerus.out
#SBATCH --error=BBmap-Phacochoerus.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=180G    # memory per cpu-core
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
...)

```bash

mkdir -p ./fastq_downloads

for accession in "${accession_list[@]}"; do
    prefetch "$accession" &
done

wait # Wait for all prefetch processes to complete

for accession in "${accession_list[@]}"; do
    # Check if the compressed fastq file already exists
    if [ ! -f "./fastq_downloads/${accession}_1.fastq.gz" ] || [ ! -f "./fastq_downloads/${accession}_2.fastq.gz" ]; then
        # If not, download and compress the fastq files
        fastq-dump --split-files --outdir ./fastq_downloads SRR26417696
        pigz ./fastq_downloads/${accession}_*.fastq
    else
        echo "Compressed files for accession ${accession} already exist, skipping."
    fi
done
```


```bash
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
srun --cpus-per-task=10 --mem=180G --time=52:00:00 --pty bash
srun --cpus-per-task=10 --mem=10G --time=24:00:00 --pty bash
find skims_processing_pipeline/ -type f -name '*.fq' -exec gzip {} \;

find /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skims_processing_pipeline/ -type f -name "*.fq" | xargs -I {} -P 30 gzip -k {} \;



cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
ln -s /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Clupea/*.fastq.gz ./
conda activate tutorial
bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testPhacochoerus -r 30 -f 30 > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testPhacochoerus.log
```

#skims_processing_pipeline.sh
#Mandatory inputs:
#-x  path to folder containing reads (split reads that need to be merged)
#-g  path to reference library
#-r  threads for RESPECT; default: 8
#-d  number of iteration cycles for RESPECT, default: 1000
#-f  number of cores for SKMER, default: 8"



```bash
srun --cpus-per-task=30 --mem=40G --time=20:00:00 --pty bash
conda activate tutorial
bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie -r 39 -f 39 > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie_11feb24_screen.log



Clupea_gskim_pocess_20nov23.sh
#!/bin/bash
#SBATCH --job-name=skim_process_clupea
#SBATCH --output=skim_process_clupea_20nov23.out
#SBATCH --error=skim_process_clupea_20nov23.err
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=18G    # memory per cpu-core
#SBATCH --time=100:00:00
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=homerejalves.monteiro@sund.ku.dk

cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
conda activate tutorial
bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea -r 56 -f 56 > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skim_process_clupea_18dec23_screen.log
```

```bash
#On a subset for Clupea dataset
#Skmer distance on existing library
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
bash ../tsv_to_phymat.sh ref-dist-mat.txt ref-dist-mat_clupea21nov.phy
# Infer backbone from scratch
../fastme-2.1.5/binaries/fastme-2.1.5-linux64  -i ref-dist-mat_clupea21nov.phy  -o backbone-fastme_clupea21nov.tre
# Recompute branch lengths:
../fastme-2.1.5/binaries/fastme-2.1.5-linux64  -i ref-dist-mat_clupea21nov.phy  -o backbone-fastme_clupea21nov_2.tre -u backbone-fastme_clupea21nov.tre
#
nw_reroot backbone-fastme_clupea21nov_2.tre |nw_display -


 /----------------------------------------+ unclassified-kra ERR9709184: Kadriorg Wreck, Estonia, 14th century
 |
=+                                        /---------------------------------------------------+ unclassified-kra ERR9709103: Giecz, Poland, 9th-13th century
 |                                        |
 \----------------------------------------+              /------------+ unclassified-kra ERR9709099: Giecz, Poland, 9th-13th century
                                          |              |
                                          |              |                             /-----------------------------------------+ unclassified-kra ERR9709108: Giecz, Poland, 9th-13th century
                                          \--------------+                         /---+
                                                         |                         |   \-----------------------------------+ unclassified-kra ERR9709094: Giecz, Poland, 9th-13th century
                                                         |                         |
                                                         \-------------------------+     /--------------------------------------------+ unclassified-kra ERR9709123: Mala Nieszawka, Poland, 14th-15th century
                                                                                   |     |
                                                                                   |     |    /----------------------------------+ unclassified-kra ERR9709092: Giecz, Poland, 9th-13th century
                                                                                   \-----+    |
                                                                                         |    |  /---------------------------------------+ unclassified-kra ERR9709199: Selso-Vestby, Denmark, 1290-1380 CE
                                                                                         \----+  |
                                                                                              |  |                                           /--------------+ unclassified-kra ERR9709180: Kolowbrzweg Budzistowo, Poland, 750-850 CE
                                                                                              |  |                                           |
                                                                                              |  /-------------------------------------------+              /--------------+ unclassified-kra ERR9709174: BAL_22_HER054, Kolowbrzweg Budzistowo, Poland, 750-850 CE
                                                                                              \--|                                           |              |
                                                                                                 |                                           \--------------+         /-------------------+ unclassified-kra ERR9709173: BAL_22_HER052b, Kolowbrzweg Budzistowo, Poland, 750-850 CE
                                                                                                 |                                                          \---------+
                                                                                                 |                                                                    \---------+ unclassified-kra ERR9709172: BAL_22_HER052a, Kolowbrzweg Budzistowo, Poland, 750-850 CE
                                                                                                 |          /----------------------------------+ unclassified-kra ERR9709155: BAL_22_HER045a,  Kolowbrzweg Budzistowo, Poland, 750-850 CE
                                                                                                 /----------+
                                                                                                 |          \-------------+ unclassified-kra ERR9709149: BAL_22_HER042a,  Kolowbrzweg Budzistowo, Poland, 750-850 CE
                                                                                                 |
                                                                                                /----------------------------------+ unclassified-kra ERR9709124: BAL_22_HER009a, Mala Nieszawka, Poland, 14th-15th century
                                                                                                ||
                                                                                                ||   /--+ unclassified-kra ERR9709109: BAL_22_HER005a,  Giecz, Poland, 9th-13th century
                                                                                                ||   |
                                                                                                |\   |        / unclassified-kra ERR9709101: BAL_22_HER003c, Giecz, Poland, 9th-13th century
                                                                                                |    |   /----+
                                                                                                |    |   |    \+ unclassified-kra ERR9709100: BAL_22_HER003b, Giecz, Poland, 9th-13th century
                                                                                                \----+   |     /--------+ unclassified-kra ERR9709107: BAL_22_HER004d, Giecz, Poland, 9th-13th century
                                                                                                     |   |  /--+
                                                                                                     |   |  |  \-----+ unclassified-kra ERR9709093: BAL_22_HER001b, Giecz, Poland, 9th-13th century
                                                                                                     |   |  |
                                                                                                     |   |  | /+ unclassified-kra ERR9709096: BAL_22_HER002b, Giecz, Poland, 9th-13th century
                                                                                                     \---+  | |
                                                                                                         |  | |      /------+ unclassified-kra ERR9709198: BAL_22_HER118b, Selso-Vestby, Denmark, 1290-1380 CE
                                                                                                         |  | |      |
                                                                                                         |  | |      |          /----------+ unclassified-kra ERR9709191: BAL_22_HER112b, Selso-Vestby, Denmark, 1290-1380 CE
                                                                                                         |  | |      | /--------+
                                                                                                         |  | |      | |        \----------------------------------------------+ unclassified-kra ERR9709190: BAL_22_HER112a, Selso-Vestby, Denmark, 1290-1380 CE
                                                                                                         |  | |     /+ |
                                                                                                         |  | |     || |               /--------+ unclassified-kra ERR9709158: BAL_22_HER045d, Selso-Vestby, Denmark, 1290-1380 CE
                                                                                                         \--+ |     || |               |
                                                                                                            | |     || |             /-/---------+ unclassified-kra ERR9709147: BAL_22_HER041a, Truso Janów Pomorski, Poland, 800-850 CE
                                                                                                            | |     || |             | |
                                                                                                            | |     |\-+             | \ /-----------+ unclassified-kra ERR9709145: BAL_22_HER039e, Truso Janów Pomorski, Poland, 800-850 CE
                                                                                                            | |     |  |             | | |
                                                                                                            | |     |  |             | \-+/--------------+ unclassified-kra ERR9709134: BAL_22_HER037b, Truso Janów Pomorski, Poland, 800-850 CE
                                                                                                            | |     |  |       /-----+   ||
                                                                                                            | |     |  |       |     |   \+    /+ unclassified-kra ERR9709132: BAL_22_HER036d, Truso Janów Pomorski, Poland, 800-850 CE
                                                                                                            | |     |  |       |     |    \----+
                                                                                                            | |     |  \-------+     |         \-------------------------+ unclassified-kra ERR9709130: BAL_22_HER036b, Truso Janów Pomorski, Poland, 800-850 CE
                                                                                                            | |     |          |     |
                                                                                                            | |     |          |     \------------------------------------+ unclassified-kra ERR9709142: BAL_22_HER039b, Truso Janów Pomorski, Poland, 800-850 CE
                                                                                                            \-+     |          |
                                                                                                              |     /          \------+ unclassified-kra ERR9709144: BAL_22_HER039d, Truso Janów Pomorski, Poland, 800-850 CE
                                                                                                              |     |
                                                                                                              |     |          /---------------------------------+ unclassified-kra ERR9709196: BAL_22_HER117, Selso-Vestby, Denmark, 1290-1380 CE
                                                                                                              |     |         /+
                                                                                                              |     |   /-----+\----------------------------------+ unclassified-kra ERR9709115: BAL_22_HER006c, Mala Nieszawka, Poland, 14th-15th century
                                                                                                              |     |   |     |
                                                                                                              |     |   | /   \----+ unclassified-kra ERR9709114: BAL_22_HER006b, Mala Nieszawka, Poland, 14th-15th century
                                                                                                              |     |   | |
                                                                                                              |     |   \-+-----+ unclassified-kra ERR9709127: BAL_22_HER009d,  Mala Nieszawka, Poland, 14th-15th century
                                                                                                              |     | /---+
                                                                                                              |     | |   |                              | unclassified-kra ERR9709136: BAL_22_HER038b, Truso Janów Pomorski, Poland, 800-850 CE
                                                                                                              |     |/+   \---------------------------------+
                                                                                                              |     |||                                     \------+ unclassified-kra ERR9709135: BAL_22_HER038a, Truso Janów Pomorski, Poland, 800-850 CE
                                                                                                              |     |||
                                                                                                              |     /+\+ unclassified-kra ERR9709125: BAL_22_HER009b, Mala Nieszawka, Poland, 14th-15th century
                                                                                                              |     ||
                                                                                                              |     \\-+ unclassified-kra ERR9709126: BAL_22_HER009c, Mala Nieszawka, Poland, 14th-15th century 
                                                                                                              |     |
                                                                                                              |     \-----+ unclassified-kra ERR9709118: BAL_22_HER007c, Mala Nieszawka, Poland, 14th-15th century  
                                                                                                              |     |
                                                                                                              |     | /-+ unclassified-kra ERR9709240: BAL_22_M-HER052, Kalmarsund, Sweden, 2010
                                                                                                              \-----+ |
                                                                                                                    | /-+ unclassified-kra ERR9709243: BAL_22_M-HER055, Kalix, Finland, 2002
                                                                                                                    | |
                                                                                                                    | /+ unclassified-kra ERR9709223: BAL_22_M-HER023, Idefjord - Inner, Norway, 2010
                                                                                                                    | \
                                                                                                                    | \ unclassified-kra ERR9709203: BAL_22_M-HER002, Karmoy, Norway, 2002
                                                                                                                    | |
                                                                                                                    | \-+ unclassified-kra ERR9709219: BAL_22_M-HER019, Kalmarsund, Sweden, 2010
                                                                                                                    | |
                                                                                                                    | \-+ unclassified-kra ERR9709218: BAL_22_M-HER018, Kalix, Finland, 2002
                                                                                                                    | |
                                                                                                                    | \-+ unclassified-kra ERR9709238: BAL_22_M-HER045, Idefjord - Outer, Norway, 2010
                                                                                                                    | |
                                                                                                                    | \-+ unclassified-kra ERR9709234: BAL_22_M-HER037, Måseskär, Sweden, 2003
                                                                                                                    | |
                                                                                                                    | \-+ unclassified-kra ERR9709235: BAL_22_M-HER038, Risor, Norway, 2003
                                                                                                                    | |
                                                                                                                    | \-+ unclassified-kra ERR9709208: BAL_22_M-HER008, More, Norway, 2002
                                                                                                                    | |
                                                                                                                    | \-+ unclassified-kra ERR9709249: BAL_22_M-HER062, Idefjord - Inner, Norway, 2010
                                                                                                                    \-+
                                                                                                                      /-+ unclassified-kra ERR9709214: BAL_22_M-HER014, Måseskär, Sweden, 2003
                                                                                                                      |
                                                                                                                      /-+ unclassified-kra ERR9709211: BAL_22_M-HER011, Risor, Norway, 2003
                                                                                                                      |
                                                                                                                    | unclassified-kra ERR9709232: BAL_22_M-HER035, Måseskär, Sweden, 2003
                                                                                                                      /
                                                                                                                      \-+ unclassified-kra ERR9709204: BAL_22_M-HER003, Karmoy, Norway, 2002
                                                                                                                      |
                                                                                                                      \-+ unclassified-kra ERR9709253: BAL_22_M-HER066, Risor, Norway, 2003

 |-----------------------------------|-----------------------------------|-----------------------------------|----------------------------------|-----------------------------------|------
 0                                0.05                                 0.1                                0.15                                0.2                                0.25
 substitutions/site



```


Magpie_gskim_pocess_20nov23.sh
#!/bin/bash
#SBATCH --job-name=skim_process_magpie
#SBATCH --output=skim_process_magpie_20nov23.out
#SBATCH --error=skim_process_magpie_20nov23.err
#SBATCH --cpus-per-task=35
#SBATCH --mem=200G    # memory per cpu-core
#SBATCH --time=300:00:00
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=homerejalves.monteiro@sund.ku.dk

cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie
conda activate tutorial
bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie -r 32 -f 32 > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/skim_process_magpie15dec23.log






# 13dec23 bwa_mem
# pwd=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd
# /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/01_bwa_mem_launch.sh

```bash
#!/bin/bash

# Clean session

rm /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/00_scripts/BWA*sh

# launch scripts for Colosse

for file in $(ls /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skims_processing_pipeline/kraken/*.fq.gz|sed -e 's/_.fq.gz//'|sort -u)

do

base=$(basename "$file")

        toEval="cat /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/01_bwa_mem.sh | sed 's/__BASE__/$base/g'"; eval $toEval > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/BWA_$base.sh
done




# For HanClupea
```bash
#!/bin/bash

# Directory containing the folders with .sra files
base_dir="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testHanClupea"

# Loop through each directory starting with SRR*
for dir in "$base_dir"/SRR*/; do
    echo "Checking directory: $dir"

    # Move into the directory, or skip if it doesn't exist
    cd "$dir" || continue

    # Check for existing .fastq.gz files in the directory
    if ls *.fastq.gz 1> /dev/null 2>&1; then
        echo "Skipped: .fastq.gz files already exist in $dir"
        # Return to the base directory and continue with next directory
        cd "$base_dir" || exit
        continue
    fi

    echo "Processing directory: $dir"

    # Enable globbing to match file names
    shopt -s nullglob
    shopt -s dotglob

    # Find .sra files in the directory
    sra_files=(./*.sra)

    if [ ${#sra_files[@]} -gt 0 ]; then
        # Convert .sra files to FASTQ and compress to .gz format
        for file in "${sra_files[@]}"; do
            accession=$(basename "$file" .sra)
            fastq-dump --split-files "$file" && pigz "${accession}"_*.fastq
        done

        # Remove .sra files after conversion
        rm -f ./*.sra
    fi

    # Return to the base directory
    cd "$base_dir" || exit
done

echo "All processing complete."
```

```bash
#!/bin/bash
#SBATCH --job-name=HanSkmin_sbatch
#SBATCH --output=Skmin_sbatch_14jan.out
#SBATCH --error=Skmin_sbatch_14jan.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=240G # memory per cpu-core
#SBATCH --time=100:00:00
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=homerejalves.monteiro@sund.ku.dk
```
Skmer's main strength lies in its ability to quickly estimate distances between lcWGS individuals, without a reference genome, reducing computing requirements and bioinformatics steps for monitoring populations. It is also a great alternative undistinguished taxa using mtDNA barcoding.

Modules (FST and MDS plots) need to incorporate statistics from Skmer to discard individuals and populations that do not meet certain thresholds for error rate, coverage, read length, genome size, and heterozygosity.

Ancient DNA specimens are not suitable for the current Skmer pipeline due to their significantly higher heterozygosity (θ).

Stronger genetic clustering was observed at sampling sites using ANGSD compared to Skmer, likely due to the removal of high LD regions.

Clear genetic breaks between ecotypes were identified by both methods. These methods highlight that many downstream applications of genome skims can be successfully performed if we have an accurate way to measure the distance between two genomes from which genome skims are generated.

Skmer can be upscaled to population diversity metrics estimation since the mean Hamming distance between samples is equivalent to heterozygosity. Fst (genetic differentiation) is defined based on the distance between samples, making Skmer well-suited for identifying structure in species at early stages of domestication and breeding.
