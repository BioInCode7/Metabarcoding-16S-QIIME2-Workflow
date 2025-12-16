#!/bin/bash

###############################################################################
# PICRUSt2 functional inference
###############################################################################
#
# IMPORTANT:
# - Run this script OUTSIDE QIIME2
# - Use a dedicated PICRUSt2 environment
# - Input table must NOT be rarefied
#
###############################################################################

set -euo pipefail

# ===============================
# INPUT FILES
# ===============================

FEATURE_TABLE_BIOM="../exported_data/picrust2/feature-table.biom"
REP_SEQS_FASTA="../exported_data/picrust2/dna-sequences.fasta"

# ===============================
# OUTPUT DIRECTORY
# ===============================

OUTPUT_DIR="./picrust2_out_pipeline"
N_THREADS=12

###############################################################################
# RUN PICRUSt2
###############################################################################

picrust2_pipeline.py \
  -s "${REP_SEQS_FASTA}" \
  -i "${FEATURE_TABLE_BIOM}" \
  -o "${OUTPUT_DIR}" \
  -p "${N_THREADS}" \
  --stratified \
  --coverage \
  --per_sequence_contrib \
  --verbose \
  --remove_intermediate

###############################################################################
# END OF PICRUSt2 PIPELINE
###############################################################################
