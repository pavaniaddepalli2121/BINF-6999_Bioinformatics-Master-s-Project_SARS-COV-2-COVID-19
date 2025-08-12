#!/usr/bin/env bash
set -euo pipefail

############################################
# SARS-CoV-2 QC Pipeline for Large Datasets
# Handles 400+ FASTQ pairs automatically
# Modified with lenient trimming params to improve retention
# Author: Pavani
# Date: $(date)
############################################

# ---- USER SETTINGS ----
DATA_DIR="/mnt/d/waste water project/filtered_sarscov2"
OUT_DIR="$DATA_DIR/qc_out"
RAW_QC_DIR="$OUT_DIR/fastqc_raw"
TRIM_QC_DIR="$OUT_DIR/fastqc_trimmed"
TRIM_DIR="$OUT_DIR/trimmed"
SUMMARY_TSV="$OUT_DIR/qc_read_counts.tsv"
THREADS=12   # adjust based on your CPU

# Path to Trimmomatic adapters (edit if needed)
ADAPTERS="$(dirname "$(which trimmomatic)")/../share/trimmomatic*/adapters/TruSeq3-PE.fa"

# ---- CREATE OUTPUT DIRECTORIES ----
mkdir -p "$RAW_QC_DIR" "$TRIM_QC_DIR" "$TRIM_DIR"

echo -e "sample\traw_reads\ttrimmed_reads\tpercent_retained" > "$SUMMARY_TSV"

count_reads() {
    local fq="$1"
    zcat -f "$fq" | awk 'END{print NR/4}'
}

# ---- LOOP OVER ALL R1 FILES ----
shopt -s nullglob
for R1 in "$DATA_DIR"/*_R1*.fastq.gz; do
    sample=$(basename "$R1" | sed -E 's/_R1.*\.fastq\.gz//')
    R2="${R1/_R1/_R2}"

    if [[ ! -f "$R2" ]]; then
        echo "WARNING: Missing R2 for $sample"
        continue
    fi

    echo ">>> Processing: $sample"

    # 1. FastQC on raw reads
    fastqc -t "$THREADS" -o "$RAW_QC_DIR" "$R1" "$R2"

    # 2. Trimming with lenient params: sliding window quality lowered to 20, min length lowered to 50
    R1P="$TRIM_DIR/${sample}_clean_R1.fastq.gz"
    R1U="$TRIM_DIR/${sample}_clean_R1.unpaired.fastq.gz"
    R2P="$TRIM_DIR/${sample}_clean_R2.fastq.gz"
    R2U="$TRIM_DIR/${sample}_clean_R2.unpaired.fastq.gz"

    trimmomatic PE -threads "$THREADS" -phred33 \
        "$R1" "$R2" \
        "$R1P" "$R1U" \
        "$R2P" "$R2U" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 SLIDINGWINDOW:4:20 MINLEN:50

    # 3. FastQC on trimmed reads
    fastqc -t "$THREADS" -o "$TRIM_QC_DIR" "$R1P" "$R2P"

    # 4. Read count summary
    raw_reads=$(count_reads "$R1")
    trimmed_reads=$(count_reads "$R1P")
    percent=$(awk -v r="$raw_reads" -v t="$trimmed_reads" 'BEGIN{if(r==0)print 0; else print (t/r)*100}')
    printf "%s\t%d\t%d\t%.2f\n" "$sample" "$raw_reads" "$trimmed_reads" "$percent" >> "$SUMMARY_TSV"
done

# 5. MultiQC summary
multiqc "$OUT_DIR" -o "$OUT_DIR/multiqc_report"

echo "QC completed for all samples!"
echo "Reports and trimmed FASTQs are in: $OUT_DIR"