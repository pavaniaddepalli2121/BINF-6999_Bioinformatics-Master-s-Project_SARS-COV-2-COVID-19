#!/usr/bin/env bash
set -euo pipefail

BAM_DIR="/mnt/d/waste water project/filtered_sarscov2"
SORTED_DIR="${BAM_DIR}/sorted_bam"
THREADS=4   # reduce threads to limit memory usage
MEM_PER_THREAD="500M"  # limit memory per thread

mkdir -p "$SORTED_DIR"

for BAM in "$BAM_DIR"/*.bam; do
  sample=$(basename "$BAM" .bam)
  sorted_bam="$SORTED_DIR/${sample}.sorted.bam"

  echo "Sorting $sample.bam ..."
  samtools sort -@ "$THREADS" -m "$MEM_PER_THREAD" -o "$sorted_bam" "$BAM"

  echo "Indexing $sample.sorted.bam ..."
  samtools index "$sorted_bam"
done

echo "All BAM files sorted and indexed in $SORTED_DIR"

