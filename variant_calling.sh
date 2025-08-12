#!/usr/bin/env bash
set -euo pipefail

REF="/mnt/d/waste water project/NC_045512.2.fasta"
SORTED_DIR="/mnt/d/waste water project/filtered_sarscov2/sorted_bam"
VCF_DIR="/mnt/d/waste water project/filtered_sarscov2/vcfs"
THREADS=4

mkdir -p "$VCF_DIR"

for BAM in "$SORTED_DIR"/*.sorted.bam; do
  sample=$(basename "$BAM" .sorted.bam)
  out_vcf="$VCF_DIR/${sample}.vcf.gz"

  echo "Calling variants for $sample..."

  bcftools mpileup -f "$REF" -a DP,AD -Q 20 -B -Ou "$BAM" | \
  bcftools call -mv -Oz -o "$out_vcf"

  bcftools index "$out_vcf"
done

echo "Variant calling completed. VCFs are in $VCF_DIR"
