# 1. Skmer-2 with spectra from these more complete genomes
```bash

module load parallel
conda activate skimming_echarvel

SANGER_GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeOedulis/ncbi_dataset/data/GCF_947568905.1/GCF_947568905.1_xbOstEdul1.1_genomic.fna"
SANGER_RESPECT_SPECTRA="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/respect-reference.txt"
PROCESS_FASTQ="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/skims_processing_pipeline_oct24/kraken_subsample/"
DATE=$(date +%d.%m)

cd /projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis

python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/skmer_new-err.py --debug reference $PROCESS_FASTQ -r $SANGER_RESPECT_SPECTRA -p 2 -l ./skmer2_sangerRespectspectra_suboedulis_$DATE
