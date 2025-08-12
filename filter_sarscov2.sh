#!/bin/bash

REF="/mnt/d/wbe_project_new/NC_045512.2.fasta"
FASTQ_DIR="/mnt/c/Users/drpav/University of Guelph/Opeyemi Lawal - Airport_WW_Sequence/airport-guelph-II"
OUT_DIR="/mnt/d/wbe_project_new/filtered_sarscov2"
mkdir -p "$OUT_DIR"

THREADS=4

for R1 in "$FASTQ_DIR"/*_R1_001.fastq.gz; do
    BASENAME=$(basename "$R1" "_R1_001.fastq.gz")
    R2="$FASTQ_DIR/${BASENAME}_R2_001.fastq.gz"

    if [ ! -f "$R2" ]; then
        echo "Skipping $BASENAME: R2 missing"
        continue
    fi

    echo "Processing $BASENAME"

    bwa mem -t $THREADS "$REF" "$R1" "$R2" | \
    samtools view -b -F 4 - | \
    samtools sort -o "$OUT_DIR/${BASENAME}_sarscov2.bam"

    samtools index "$OUT_DIR/${BASENAME}_sarscov2.bam"

    samtools fastq "$OUT_DIR/${BASENAME}_sarscov2.bam" \
      -1 "$OUT_DIR/${BASENAME}_sarscov2_R1.fastq.gz" \
      -2 "$OUT_DIR/${BASENAME}_sarscov2_R2.fastq.gz" \
      -0 /dev/null -s /dev/null -n

    echo "Finished $BASENAME"
done