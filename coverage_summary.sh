#!/bin/bash

BAM_DIR="/mnt/d/wbe_project_new/filtered_sarscov2"
GENOME_LENGTH=29903
OUT_CSV="/mnt/d/wbe_project_new/coverage_summary_extended.csv"

# Write CSV header
echo "Sample,Total_Mapped_Reads,Coverage_1x_Percentage,Coverage_10x_Percentage" > "$OUT_CSV"

for BAM in "$BAM_DIR"/*_sarscov2.bam; do
    SAMPLE=$(basename "$BAM" "_sarscov2.bam")
    echo "Processing sample: $SAMPLE"

    TOTAL_MAPPED=$(samtools idxstats "$BAM" | awk '{mapped+=$3} END {print mapped}')

    # Coverage positions with depth >=1
    POS_1X=$(samtools depth "$BAM" | awk '{if($3>=1) count++} END {print count}')
    POS_1X=${POS_1X:-0}

    # Coverage positions with depth >=10
    POS_10X=$(samtools depth "$BAM" | awk '{if($3>=10) count++} END {print count}')
    POS_10X=${POS_10X:-0}

    COVERAGE_1X=$(awk -v pos="$POS_1X" -v len="$GENOME_LENGTH" 'BEGIN {printf "%.2f", (pos/len)*100}')
    COVERAGE_10X=$(awk -v pos="$POS_10X" -v len="$GENOME_LENGTH" 'BEGIN {printf "%.2f", (pos/len)*100}')

    echo "$SAMPLE,$TOTAL_MAPPED,$COVERAGE_1X,$COVERAGE_10X" >> "$OUT_CSV"
done

echo "Extended coverage summary saved to $OUT_CSV"
