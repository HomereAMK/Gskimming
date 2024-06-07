#!/bin/bash

# Clean session

rm /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/00_scripts/FULLANGSD*sh

# launch scripts for Colosse

for file in $(ls /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/subsampled_4/*.fastq |sed -e 's/__merged.fastq//'| sed -e 's/unclassified-kra_unclassified-seq_//'| sort -u)

do

base=$(basename "$file")

        toEval="cat /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/4x_fullAngsdProcess.sh | sed 's/__BASE__/$base/g'"; eval $toEval > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/4x_FULLANGSD_$base.sh
done


#Submit jobs
for i in $(ls /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/scripts/4x_FULLANGSD*sh); do sbatch $i; done
