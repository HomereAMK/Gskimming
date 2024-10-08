
```bash
INPUT="/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/mappability/fastq"

        bwa mem -t 30 -R "@RG\tID:$base\tSM:$base\tPL:Illumina" "$GENOME" "$filepath1" "$filepath2" > "$DATAOUTPUT/$base.sam" 2> "$DATAOUTPUT/$base.bwa.log"




# istats 

scale is 1000
Number of query intervals (n) is 116
Number of reference intervals (m) is 101
Overlap is 54
eta is 0.00102011
DP p-value is 8.66247e-06
PB p-value is 2.63722e-10