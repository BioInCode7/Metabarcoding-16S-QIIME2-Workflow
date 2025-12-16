# NOTE:
# Truncation lengths were chosen based on empirical ASV length distributions
# obtained from DADA2 denoising (mean ~407 bp). Multiple classifiers were trained
# to assess robustness and ensure length-matching with inferred ASVs.


#!/bin/bash
set -euo pipefail

###############################################################################
# Train SILVA 138 Naive Bayes classifier for 16S V3–V4 (length-matched to ASVs)
###############################################################################
#
# This classifier is trained to match:
#   - Amplicon region: V3–V4
#   - Primer pair: 341F / 806R
#   - ASV length inferred by DADA2 (~407 bp)
#
# Final truncation length selected: 400 bp
#
###############################################################################

# ===============================
# BASE DIRECTORIES
# ===============================

REPO_BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

SOURCE_DIR="${REPO_BASE_DIR}/02_qiime2_pipeline/classifiers/silva-138-source"
OUT_DIR="${REPO_BASE_DIR}/02_qiime2_pipeline/classifiers/silva-138-v3v4-420bp"

SEQS_QZA="${SOURCE_DIR}/silva-138-99-seqs.qza"
TAX_QZA="${SOURCE_DIR}/silva-138-99-tax.qza"

# ===============================
# Primers (V3–V4)
# ===============================

FORWARD_PRIMER="CCTACGGGNGGCWGCAG"
REVERSE_PRIMER="GACTACHVGGGTATCTAATCC"

TRUNC_LEN=420
N_JOBS=16

mkdir -p "${OUT_DIR}"

echo "============================================================"
echo "Training SILVA 138 V3–V4 classifier (${TRUNC_LEN} bp)"
echo "============================================================"

# ===============================
# Step 1 — Extract reads
# ===============================

qiime feature-classifier extract-reads \
  --i-sequences "${SEQS_QZA}" \
  --p-f-primer "${FORWARD_PRIMER}" \
  --p-r-primer "${REVERSE_PRIMER}" \
  --p-trunc-len "${TRUNC_LEN}" \
  --p-n-jobs "${N_JOBS}" \
  --o-reads "${OUT_DIR}/silva-138-V3V4-420-ref-seqs.qza"

# ===============================
# Step 2 — Train classifier
# ===============================

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads "${OUT_DIR}/silva-138-V3V4-420-ref-seqs.qza" \
  --i-reference-taxonomy "${TAX_QZA}" \
  --o-classifier "${OUT_DIR}/silva-138-V3V4-420-classifier.qza"

echo "✅ Classifier trained successfully:"
echo "   ${OUT_DIR}/silva-138-V3V4-420-classifier.qza"