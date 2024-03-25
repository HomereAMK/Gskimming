# Gskimming
 Analyses of Nuclear Reads Obtained using Genome skimming

# Module list in Mjolnir 
```  1) skmer/3.3.0   2) apples/20231016   3) fasttree/2.1.11   4) mash/2.3   5) fastme/2.1.6.1  /projects/mjolnir1/people/sjr729/softwares/guppy-64 ```

#if CONSULT used, it needs a lot of memory
```
#salloc --qos=debug --nodes=1 -c 20 --mem-per-cpu 18000 -t 140000
srun --cpus-per-task=40 --mem=200g --time=300:00:00 --pty bash
srun --cpus-per-task=30 --mem=70g --time=48:00:00 --pty bash


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
        fasterq-dump --split-files --outdir ./fastq_downloads "$accession"
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



##### Clupea redl fastq
```bash
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testHanClupea
#!/bin/bash

# Loop through each directory starting with SRR*
for dir in SRR*/; do
    cd "$dir" || continue  # Move into the directory, or skip if it doesn't exist

    # Check if there are .sra files in the directory
    shopt -s nullglob
    shopt -s dotglob

    sra_files=(SRR*.sra)  # Find .sra files

    if [ ${#sra_files[@]} -gt 0 ]; then
        # Convert .sra files to FASTQ
        for file in *.sra; do
            fastq-dump --split-files "$file"
        done

        # Remove .sra files
        rm -f *.sra
        pigz ./fastq_downloads/${accession}_*.fastq

    fi

    # Check if the directory is empty after removing .sra files
    if [ ! "$(ls -A .)" ]; then
        cd ..
        rm -r "$dir"  # Delete the directory if it's empty
    else
        cd ..
    fi
done
```
```bash
conda activate tutorial
bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea -r 32 -f 32 > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skim_process_clupea_11dec23_screen.log


srun --cpus-per-task=60 --mem=200G --time=72:00:00 --pty bash
conda activate tutorial
bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea -r 56 -f 56 > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skim_process_clupea_13dec23_screen.log
```


bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie -r 39 -f 39 > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testMagpie/skim_process_magpie_09jan24_screen.log









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


#Submit jobs
for i in $(ls /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/BWA*sh); do sbatch $i; done
```


# /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/01_bwa_mem.sh
```bash
#!/bin/bash
#SBATCH --job-name=__BASE__bwa_clupea
#SBATCH --output=98_log_files/Map__BASE__bwa_clupea.out
#SBATCH --error=98_log_files/Map__BASE__bwa_clupea.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=25g   
#SBATCH --time=15:00:00
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=homerejalves.monteiro@sund.ku.dk

# Load all required modules for the job
module load gcc/8.2.0
module load bwa/0.7.17
module load samtools/1.18

# Global variables
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/mapped"
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skims_processing_pipeline/kraken/"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
base=__BASE__


#     Align reads
    echo "Aligning $base"
    ID=$(echo "@RG\tID:$base\tSM:$base\tPL:Illumina")

  # Align reads 1 step
    bwa mem -t "$NCPU" \
        -R "$ID" \
        "$GENOME" \
        "$DATAINPUT"/"$base"_fq.gz >"$DATAOUTPUT"/"$base".sam


        # Create bam file
    echo "Creating bam for $base"

    samtools view -bS -h -q 20 -F 4 \
    "$DATAOUTPUT"/"$base".sam >"$DATAOUTPUT"/"$base".bam


     echo "Creating sorted bam for $base"
        samtools sort "$DATAOUTPUT"/"$base".bam -o "$DATAOUTPUT"/"$base".sort.minq20.bam
        samtools index "$DATAOUTPUT"/"$base".sort.minq20.bam
   
   # Clean up
    echo "Removing "$DATAOUTPUT"/"$base".sam"
    echo "Removing "$DATAOUTPUT"/"$base".bam"

        rm "$DATAOUTPUT"/"$base".sam
        rm "$DATAOUTPUT"/"$base".bam

```
# Load all required modules for the job
module load gcc/8.2.0
module load bwa/0.7.17
module load samtools/1.18

# Global variables
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/mapped"
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skims_processing_pipeline/kraken/"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"


#     Align reads
    echo "Aligning $base"
    ID=$(echo "@RG\tID:$base\tSM:$base\tPL:Illumina")

  # Align reads 1 step
    bwa mem -t "$NCPU" \
        -R "$ID" \
        "$GENOME" \
        "$DATAINPUT"/"$base"_fq.gz >"$DATAOUTPUT"/"$base".sam


        # Create bam file
    echo "Creating bam for $base"

    samtools view -bS -h -q 20 -F 4 \
    "$DATAOUTPUT"/"$base".sam >"$DATAOUTPUT"/"$base".bam


     echo "Creating sorted bam for $base"
        samtools sort "$DATAOUTPUT"/"$base".bam -o "$DATAOUTPUT"/"$base".sort.minq20.bam
        samtools index "$DATAOUTPUT"/"$base".sort.minq20.bam

   # Clean up
    echo "Removing "$DATAOUTPUT"/"$base".sam"
    echo "Removing "$DATAOUTPUT"/"$base".bam"

        rm "$DATAOUTPUT"/"$base".sam
        rm "$DATAOUTPUT"/"$base".bam


# 18dec MarkDuplicates launcher 02_MarkDuplicates_clipOverlap_launch.sh
```bash
#!/bin/bash

# Clean session

rm 00_scripts/MARKDUP*sh

# launch scripts for Colosse

for file in $(ls /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/mapped/*sort*.bam|sed -e 's/.sort.minq20.bam//g'|sort -u)
do

base=$(basename "$file")

	toEval="cat /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/02_MarkDuplicates_clipOverlap.sh | sed 's/__BASE__/$base/g'"; eval $toEval > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/MARKDUP_$base.sh
done


#Submit jobs
for i in $(ls /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/MARKDUP*sh); do sbatch $i; done
```

# 18dec MarkDuplicates +realigned code 02_MarkDuplicates_clipOverlap.sh
```bash
#!/bin/bash
#SBATCH --job-name=markdupnClip__BASE__
#SBATCH --output=98_log_files/DuplicateClip/__BASE__markdnClip.out
#SBATCH --error=98_log_files/DuplicateClip/__BASE__markdnClip.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=25g   
#SBATCH --time=15:00:00
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=homerejalves.monteiro@sund.ku.dk


# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/

### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)


#loading modules
module load jdk/1.8.0_291  picard/2.27.5  parallel/20230822   java-jdk/8.0.112   bamutil/1.0.15  python/3.11.3   openjdk/17.0.8  samtools/1.18
module load  jdk/1.8.0_291   picard/2.27.5  parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15 python/3.11.3 openjdk/17.0.8 samtools/1.18 gatk/3.8

# Global variables
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/dedup"
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/mapped"
base=__BASE__

#tryout with NO CIGAR on MarkDuplicates
picard MarkDuplicates \
I="$DATAINPUT"/"$base".sort.minq20.bam \
O="$DATAOUTPUT"/"$base".nocig.dedup.minq20.bam \
M="$DATAOUTPUT"/"$base".duprmmetrics.txt \
REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT

wait 

#scripts ClipOverlap with NO CIGAR on MarkDuplicates
/projects/mjolnir1/apps/conda/bamutil-1.0.15/bin/bam clipOverlap \
--in "$DATAOUTPUT"/"$base".nocig.dedup.minq20.bam \
--out "$DATAOUTPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam \
--stats

wait

#ressources
module load samtools/1.12
module load parallel/20160822
module load java/1.8.0
module load bamutil/1.0.14
module load gatk/3.8-0
module load jre/1.8.0-openjdk
module load picard-tools/2.25.2

# Global variables
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/realigned"
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/dedup"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"

#move to present working dir
cd $PBS_O_WORKDIR

base=__BASE__

#scripts
#fasta seq dictionary file ref picard
#picard CreateSequenceDictionary R= $GENOME
#-fai ref
#samtools faidx $GENOME

# Index bam files
samtools index "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam 

## Create list of potential in-dels nocig
java -jar /projects/mjolnir1/apps/conda/gatk-3.8/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $GENOME \
-I "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam  \
-o "$DATAOUTPUT"/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals 

wait

## Run the indel realigner tool nocig
java -jar /projects/mjolnir1/apps/conda/gatk-3.8/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $GENOME \
-I "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam \
-targetIntervals "$DATAOUTPUT"/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals \
--consensusDeterminationModel USE_READS --nWayOut _minq20.nocig.realigned.bam

##
mv *realigned.bam realigned/
mv *realigned.bai realigned/

```















# Full launcher
```bash
#!/bin/bash

# Clean session

rm /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/scripts/FULLANGSDPROCESS*sh

# launch scripts for Colosse

for file in $(ls /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skims_processing_pipeline/consult/*__merged|sed -e 's/__merged//'|sort -u)

do

base=$(basename "$file")

        toEval="cat /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/full_AngsdProcess.sh | sed 's/__BASE__/$base/g'"; eval $toEval > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/FULLANGSDPROCESS_$base.sh
done


#Submit jobs
for i in $(ls /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/FULLANGSDPROCESS*sh); do sbatch $i; done
```

# Full angsd preprocess
```bash
#!/bin/bash
#SBATCH --job-name=__BASE__FULLANGSDPROCESS
#SBATCH --output=98_log_files/__BASE__fullangsdprocees.out
#SBATCH --error=98_log_files/__BASE__fullangsdprocees.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=25g   
#SBATCH --time=15:00:00
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=homerejalves.monteiro@sund.ku.dk

#pwd
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/
#fasta seq dictionary file ref picard
#picard CreateSequenceDictionary R= $GENOME
#-fai ref
#samtools faidx $GENOME


#module
module load gcc/8.2.0 bwa/0.7.17 samtools/1.18  jdk/1.8.0_291  picard/2.27.5  parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15 python/3.11.3 openjdk/17.0.8 gatk/3.8 

module load gcc/8.2.0 bwa/0.7.17 samtools/1.18  jdk/1.8.0_291  picard/2.27.5  parallel/20230822 java-jdk/8.0.112 bamutil/1.0.15 python/3.11.3 openjdk/17.0.8  gatk-framework/


# Global variables
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/mapped"
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skims_processing_pipeline/consult/"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
base=__BASE__

### BWA ###
#     Align reads
    echo "Aligning $base"
    ID=$(echo "@RG\tID:$base\tSM:$base\tPL:Illumina")

  # Align reads 1 step
    bwa mem -t "$NCPU" \
        -R "$ID" \
        "$GENOME" \
        "$DATAINPUT"/"$base"__merged > "$DATAOUTPUT"/"$base".sam


        # Create bam file
    echo "Creating bam for $base"

    samtools view -bS -h -q 20 -F 4 \
    "$DATAOUTPUT"/"$base".sam >"$DATAOUTPUT"/"$base".bam


     echo "Creating sorted bam for $base"
        samtools sort "$DATAOUTPUT"/"$base".bam -o "$DATAOUTPUT"/"$base".sort.minq20.bam
        samtools index "$DATAOUTPUT"/"$base".sort.minq20.bam
   
   # Clean up
    echo "Removing "$DATAOUTPUT"/"$base".sam"
    echo "Removing "$DATAOUTPUT"/"$base".bam"

        rm "$DATAOUTPUT"/"$base".sam
        rm "$DATAOUTPUT"/"$base".bam

wait

### Mark 'n' clip duplicates ###
# Global variables
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/dedup"
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/mapped"

#tryout with NO CIGAR on MarkDuplicates
picard MarkDuplicates \
I="$DATAINPUT"/"$base".sort.minq20.bam \
O="$DATAOUTPUT"/"$base".nocig.dedup.minq20.bam \
M="$DATAOUTPUT"/"$base".duprmmetrics.txt \
REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT

wait 

#scripts ClipOverlap with NO CIGAR on MarkDuplicates
/projects/mjolnir1/apps/conda/bamutil-1.0.15/bin/bam clipOverlap \
--in "$DATAOUTPUT"/"$base".nocig.dedup.minq20.bam \
--out "$DATAOUTPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam \
--stats

wait


### Realign around indels ###
# Global variables
DATAOUTPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/realigned"
DATAINPUT="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/dedup"
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"

#move to present working dir
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/


# Index bam files
samtools index "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam 

wait 

## Create list of potential in-dels nocig
java -jar /projects/mjolnir1/apps/conda/gatk-3.8/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $GENOME \
-I "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam  \
-o "$DATAOUTPUT"/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals 

wait

## Run the indel realigner tool nocig
java -jar /projects/mjolnir1/apps/conda/gatk-3.8/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $GENOME \
-I "$DATAINPUT"/"$base".nocig.dedup_clipoverlap.minq20.bam \
-targetIntervals "$DATAOUTPUT"/"$base".all_samples_for_indel_realigner.nocig.minq20.intervals \
--consensusDeterminationModel USE_READS --nWayOut _minq20.nocig.realigned.bam

##
mv *realigned.bam realigned/
mv *realigned.bai realigned/


```




```bash
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
conda activate tutorial
bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/fastq_downloads -r 36 -f 56 > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/Sprocess_clupeaAtmore_14djan23_screen.log
```


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


cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testHanClupea
conda activate tutorial
bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testHanClupea/fastq_downloads  -r 38 -f 38 > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testHanClupea/Skmin_sbatch_16jan.log

sbatch --job-name=4xSkmin_sbatch --output=4xSkmin_sbatch_16jan.out --error=4xSkmin_sbatch_16jan.err --ntasks=1 --cpus-per-task=40 --mem=180G --time=100:00:00 --mail-type=begin --mail-type=end --mail-type=fail --mail-user=homerejalves.monteiro@sund.ku.dk --wrap="cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea && conda activate tutorial && bash ../skims_processing_pipeline.sh -x /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea -r 40 -f 40 "
