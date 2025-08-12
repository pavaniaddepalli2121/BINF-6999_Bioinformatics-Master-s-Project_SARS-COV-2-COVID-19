#!/bin/bash

# Paths
BAM_DIR="/mnt/d/waste water project/filtered_sarscov2/sorted_bam"
REF="/mnt/d/waste water project/NC_045512.2.fasta"
CONS_DIR="/mnt/d/waste water project/consensus"

# Create output directory if it doesn't exist
mkdir -p "$CONS_DIR"

echo "Generating consensus FASTA sequences..."

# Get list of BAM files
BAM_FILES=("$BAM_DIR"/*.sorted.bam)

for BAM in "${BAM_FILES[@]}"; do
    sample=$(basename "$BAM" .sorted.bam)

    if [[ ! -f "$BAM" ]]; then
        echo "BAM file $BAM does not exist, skipping..."
        continue
    fi

    echo "Generating consensus for sample $sample"

    samtools mpileup -aa -A -d 1000000 -B -Q 0 -f "$REF" "$BAM" | \
        ivar consensus -p "$CONS_DIR/${sample}_consensus" -q 20 -t 0.6 -m 10
done

echo "Consensus generation complete."