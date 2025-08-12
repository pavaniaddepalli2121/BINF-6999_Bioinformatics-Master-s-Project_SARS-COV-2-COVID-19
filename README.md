# THE EVOLUTION OF SARS-COV-2 (COVID-19) VARIANTS: ANALYSIS OF LARGE GENOMIC DATASETS

**Author:** Pavani Addepalli  
**Advisors:**  Dr. Ryan Gregory and Dr. Gurjit Randhawa
**Program:** Master of Bioinformatics, University of Guelph  
**Contact:** paddepal@uoguelph.ca  

---

## 1. Project Overview
This project analyzes SARS-CoV-2 genomic data from wastewater samples collected at Toronto Pearson International Airport between November 2022 and October 2023.  
The objective is to detect cryptic and novel mutations, track variant introductions, and compare patterns to previous studies (Gregory et al., 2022; Suarez et al., 2025).  

Main goals:
- Process and quality check sequencing data.
- Perform variant calling and filtering.
- Assign lineages using Pangolin and Nextclade.
- Analyze genomic diversity with Shannon entropy.
- Identify cryptic mutations using a mutation panel.
---
## 2. Directory Structure
 results/
│ ├── qc/ # QC reports (FastQC, MultiQC)
│ ├── align/ # BAM and index files
│ ├── variants/ # VCF files (raw and filtered)
│ ├── consensus/ # Consensus genome FASTAs
│ ├── pangolin/ # Lineage assignment results
│ ├── nextclade/ # QC and mutation outputs
│ ├── phylo/ # Multiple sequence alignments and trees
│ └── figures/ # Plots and heatmaps
├── scripts/ # Pipeline scripts

---
## 3. Requirements
- **Operating System:** Linux, macOS, or WSL on Windows
- **Software:**
  - FastQC ≥ 0.11
  - Trimmomatic ≥ 0.39 (or fastp)
  - MultiQC ≥ 1.14
  - BWA ≥ 0.7.17
  - samtools ≥ 1.17
  - bcftools ≥ 1.17
  - Pangolin ≥ 4.x
  - Nextclade CLI ≥ 3.x
  - MAFFT ≥ 7.5
  - IQ-TREE2 ≥ 2.2
  - R ≥ 4.2 (tidyverse, data.table, Biostrings, seqinr)

Install all dependencies using:
```bash
mamba env create -f env/environment.yml
conda activate wastewater-sarscov2
---
## 4. Main Analysis Steps

Quality Control: FastQC and MultiQC reports.

Alignment: Map reads to NC_045512.2 reference genome using BWA.

Variant Calling: bcftools mpileup and call.

Filtering: Apply depth, quality, and allele frequency thresholds.

Consensus Generation: Create per-sample consensus FASTA.

Lineage Assignment: Run Pangolin and Nextclade.

Phylogenetic Analysis: Multiple sequence alignment with MAFFT, tree building with IQ-TREE2.

Entropy Analysis: Calculate Shannon entropy in R.

Cryptic Mutation Detection: Screen using Suarez et al. (2025) panel.
---
## 5.Key Outputs
Mutation Heatmaps – Frequency and distribution of cryptic sites.

Entropy Plots – Genome-wide diversity.

Lineage Tables – Pangolin and Nextclade results.

Phylogenetic Tree – Variant relationships.
---
## 6. References

Gregory et al. (2022) – Wastewater tracking of cryptic SARS-CoV-2 lineages.

Suarez et al. (2025) – Global cryptic mutation panel analysis.

NCBI Reference Genome NC_045512.2 – Wuhan-Hu-1.
