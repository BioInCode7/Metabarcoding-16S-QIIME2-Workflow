#!/bin/bash
set -euo pipefail

###############################################################################
# NOTE ON DATA REDUCTION
#
# This script performs random subsampling of paired-end FASTQ files
# to a fixed number of reads per sample.
#
# PURPOSE:
# - Speed up pipeline execution
# - Enable reproducible example analyses
# - Ensure equal sequencing depth across samples
#
# IMPORTANT:
# Subsampled data are intended for demonstration and testing ONLY.
# Full datasets should be used for final biological interpretation.
#
###############################################################################

# =========================
# CONFIGURATION
# =========================

INPUT_DIR="/home/jesus/Bioinformatica/Bioinfo_BIO175/JesusInsectos/Insectos_seq_2025/Gm/Gm_16S"
OUTPUT_DIR="/home/jesus/Bioinformatica/Bioinfo_BIO175/GitHUbRepoMetabarcoding/01_input_example/fastq_subsampled"

SUBSAMPLE_N=20000
SEED=42

mkdir -p "${OUTPUT_DIR}"

echo "Subsampling to ${SUBSAMPLE_N} paired-end reads per sample"
echo "Input:  ${INPUT_DIR}"
echo "Output: ${OUTPUT_DIR}"
echo "Seed:   ${SEED}"

# =========================
# SUBSAMPLING
# =========================

for R1 in "${INPUT_DIR}"/*_R1_001.fastq.gz; do
  BASENAME=$(basename "${R1}" _R1_001.fastq.gz)
  R2="${INPUT_DIR}/${BASENAME}_R2_001.fastq.gz"

  echo "▶ ${BASENAME}"

  seqtk sample -s${SEED} "${R1}" ${SUBSAMPLE_N} | gzip > \
    "${OUTPUT_DIR}/${BASENAME}_R1_sub.fastq.gz"

  seqtk sample -s${SEED} "${R2}" ${SUBSAMPLE_N} | gzip > \
    "${OUTPUT_DIR}/${BASENAME}_R2_sub.fastq.gz"
done

echo "✅ Subsampling completed"
