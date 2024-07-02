## Converting BAM to FASTQ with Bedtools and GNU Parallel


# Set locale settings
export LANGUAGE=en_DK.UTF-8
export LC_ALL=en_DK.UTF-8
export LC_CTYPE=en_DK.UTF-8
export LANG=en_DK.UTF-8

# module Bedtools and  Parallel
module load bedtools parallel

# code
```bash
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd #change pather
# Create the output directory if it doesn't exist
output_dir="mapped_fq"
mkdir -p $output_dir

# Change to the directory containing the .bam files
cd mapped

# Define the conversion function with detailed logging
convert_bam_to_fastq() {
    local bam_file=$1
    fq_file="${bam_file%.bam}.fq"
    output_path="../$output_dir/$fq_file"
    start_time=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$start_time] Starting conversion of $bam_file to $output_path" >> ../$output_dir/conversion.log
    bedtools bamtofastq -i $bam_file -fq $output_path
    end_time=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$end_time] Finished conversion of $bam_file to $output_path" >> ../$output_dir/conversion.log
}

# Loop 
for bam_file in unclassified-seq_*.bam; do
    fq_file="${bam_file%.bam}.fq"
    output_path="../$output_dir/$fq_file"
    start_time=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$start_time] Starting conversion of $bam_file to $output_path" >> ../$output_dir/conversion.log
    bedtools bamtofastq -i $bam_file -fq $output_path
    end_time=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$end_time] Finished conversion of $bam_file to $output_path" >> ../$output_dir/conversion.log
done

export -f convert_bam_to_fastq

# Find all BAM files and process them in parallel
find . -name 'unclassified-seq_*.bam' | parallel convert_bam_to_fastq

```


## Use run_skmer.sh to get distance matrix of skimmed preprocessed read mapped to the ref genome
```bash
conda activate skmer_2_test
cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea
module load parallel
DATE=$(date +%d.%m)
 /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/run_skmer.sh -i /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/angsd/mapped_fq -t 19 -p 10 -o library_runskmer_mappedreads_skmer1_condaskmer2test_$DATE
 # result are in library_runskmer_mappedreads_skmer1_condaskmer2test_26.06
 ```
