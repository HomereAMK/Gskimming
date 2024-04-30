#!/bin/bash
#SBATCH --job-name=ClupeaEstimate_0.5x
#SBATCH --output=Estimate_0.5x_%j.out
#SBATCH --error=Estimate_0.5x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=180G
#SBATCH --time=60:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=homerejalves.monteiro@sund.ku.dk

# Properly initialize Conda and activate the environment
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/ 
conda activate Mar_skmer_pip

# Load necessary modules
module load parallel kraken2

# Define DATE
DATE=$(date +%d.%m)

# Execute the script 2x 
bash subsample_and_estimate.sh -i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skims_processing_pipeline_jan24/kraken -c 0.5 -t 40 -o /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/OUT_subsample_and_estimate_0.5x > /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/${DATE}_testClupea_subsample_sbatch_0.5x.log

