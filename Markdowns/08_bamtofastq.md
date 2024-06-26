## Converting BAM to FASTQ with Bedtools and GNU Parallel



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

# Function to convert a single BAM file to FASTQ
convert_bam_to_fastq() {
    local bam_file=$1
    fq_file="${bam_file%.bam}.fq"
    output_path="../$output_dir/$fq_file"
    bedtools bamtofastq -i $bam_file -fq $output_path
}

export -f convert_bam_to_fastq

# Find all BAM files and process them in parallel
find . -name 'unclassified-seq_*.bam' | parallel convert_bam_to_fastq

```