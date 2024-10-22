#!/usr/bin/env bash

# $1 mate1
# $2 mate2
# $3 output

if [ $# -lt 3 ]; then
    echo "USAGE: $0 [first fastq(.gz) input] [second fastq(.gz) input] [output filename] ([TMPDIR])"
    exit 1
fi

set -e
set -x

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

bbmapdir=$SCRIPT_DIR
dukmem=16
splitsize=5000000  # Adjusted for testing; you can increase this based on your system's capacity

if [ $# -gt 3 ]; then
    export TMPDIR=$4
else
    export TMPDIR=.
fi

mate_1=$1
mate_2=$2

split=$(mktemp -d)

# Adjust splitsize for FASTQ records (4 lines per record)
zcat -f ${mate_1} | split -l $(( ${splitsize} * 4 )) -d - ${split}/split_R1_
zcat -f ${mate_2} | split -l $(( ${splitsize} * 4 )) -d - ${split}/split_R2_

for x in `ls ${split}/split_R1_*`; do
    # Adapter removal
    $bbmapdir/bbmap/bbduk.sh -Xmx${dukmem}g -Xms${dukmem}g \
        in1=$x in2=${x/_R1_/_R2_} \
        out1=${x/_R1_/_R1DUK_} out2=${x/_R1_/_R2DUK_} \
        ref=adapters,phix ktrim=r k=23 mink=11 hdist=1 tpe tbo overwrite=true

    rm ${x/_R1_/_R2_} $x

    # Interleave the adapter-trimmed reads
    $bbmapdir/bbmap/reformat.sh \
        in1=${x/_R1_/_R1DUK_} in2=${x/_R1_/_R2DUK_} \
        out=${x/_R1_/_INTERLEAVED_} overwrite=true

    rm ${x/_R1_/_R1DUK_} ${x/_R1_/_R2DUK_}

    # Deduplication using dedupe.sh on interleaved reads
    $bbmapdir/bbmap/dedupe.sh -Xmx${dukmem}g -Xms${dukmem}g \
        in=${x/_R1_/_INTERLEAVED_} \
        out=${x/_R1_/_DEDUP_INTERLEAVED_} \
        overwrite=true

    rm ${x/_R1_/_INTERLEAVED_}

    # **De-interleave** the deduplicated reads back into R1 and R2
    $bbmapdir/bbmap/reformat.sh \
        in=${x/_R1_/_DEDUP_INTERLEAVED_} \
        out1=${x/_R1_/_R1DEDUP_} out2=${x/_R1_/_R2DEDUP_} overwrite=true

    rm ${x/_R1_/_DEDUP_INTERLEAVED_}

    # **Interleave** the deduplicated reads properly
    $bbmapdir/bbmap/reformat.sh \
        in1=${x/_R1_/_R1DEDUP_} in2=${x/_R1_/_R2DEDUP_} \
        out=${x/_R1_/_DEDUP_CORRECTLY_INTERLEAVED_} overwrite=true

    rm ${x/_R1_/_R1DEDUP_} ${x/_R1_/_R2DEDUP_}

done

# Concatenate all the properly interleaved deduplicated files to create a single output file
cat ${split}/split_DEDUP_CORRECTLY_INTERLEAVED_* > $3

# Clean up temporary files
rm -rf "$split"
