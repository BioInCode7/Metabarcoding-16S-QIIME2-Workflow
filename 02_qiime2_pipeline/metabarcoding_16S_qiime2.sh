#!/bin/bash

###############################################################################
# Metabarcoding 16S Workflow — QIIME2
###############################################################################
#
# AUTHOR:
#   Jesús Salinas (University of Almería)
#
# REPOSITORY:
#   https://github.com/BioInCode7/Metabarcoding-16S-QIIME2-Workflow
#
# DESCRIPTION:
#   This script implements a complete 16S rRNA metabarcoding pipeline using
#   QIIME2, starting from paired-end FASTQ files and ending with:
#     - ASV inference (DADA2)
#     - Taxonomic assignment
#     - Alpha and beta diversity analyses
#     - Differential abundance testing
#     - Export of results for downstream analysis in R
#
#   The workflow is designed to be:
#     - Reproducible
#     - Transparent
#     - Conservative in interpretation
#
#   This script is intentionally verbose and heavily commented.
#   It is meant to be READ as much as it is meant to be RUN.
#
###############################################################################
#
# ⚠️ IMPORTANT CONCEPTUAL NOTES
#
# 1. This pipeline analyzes 16S rRNA amplicon data.
#    It provides information about:
#       - Community composition
#       - Relative abundance
#       - Diversity patterns
#
#    It DOES NOT measure:
#       - Metabolic activity
#       - Gene expression
#       - Enzymatic function
#
# 2. Functional inference tools (e.g. PICRUSt2) predict POTENTIAL functions
#    based on reference genomes.
#
#    → Predicted function ≠ measured function
#
# 3. All statistical outputs must be interpreted in light of:
#       - Experimental design
#       - Number of biological replicates
#       - Sequencing depth
#
###############################################################################
#
# ⚙️ WORKFLOW OVERVIEW
#
#   0. Environment setup and input preparation
#   1. Adapter and primer trimming (cutadapt)
#   2. Paired-end merging (FLASH)
#   3. Quality control and denoising (DADA2)
#   4. Feature table filtering
#   5. Taxonomic classification (Naive Bayes)
#   6. Phylogenetic tree construction
#   7. Alpha diversity analysis
#   8. Beta diversity analysis and PERMANOVA
#   9. Differential abundance analysis
#  10. Data export for R and PICRUSt2
#
###############################################################################
#
# 📌 HOW TO USE THIS SCRIPT
#
# - This script is designed so that ONLY input paths need to be modified.
# - All parameters are explicitly declared and justified.
# - Users are strongly encouraged to:
#     - Inspect intermediate outputs (.qzv)
#     - Adjust parameters based on data quality
#
# ❗ This script should NOT be treated as a black box.
#
###############################################################################
#
# 📁 EXPECTED INPUTS
#
# - Paired-end FASTQ files (Phred33)
# - Sample metadata file (QIIME2-compatible TSV)
# - Pre-trained or custom-trained taxonomic classifier
#
###############################################################################
#
# 📦 SOFTWARE REQUIREMENTS
#
# - QIIME2 (tested with versions ≥ 2022)
# - FLASH
# - biom-format
# - Conda environment properly configured
#
###############################################################################
#
# 🧠 PHILOSOPHY OF THIS PIPELINE
#
# This workflow prioritizes:
#   - Interpretability over automation
#   - Explicit decisions over defaults
#   - Methodological clarity over speed
#
# It is intentionally conservative.
#
###############################################################################

set -euo pipefail

###############################################################################
# START OF PIPELINE
###############################################################################

###############################################################################
# CONFIGURATION: PROJECT STRUCTURE AND INPUT FILES
###############################################################################
#
# This section defines all project-specific paths.
# Users should ONLY modify this block to adapt the pipeline to their system.
#
# Best practice:
# - Keep raw data, intermediate files and results clearly separated
# - Avoid relative paths inside the pipeline logic
#
###############################################################################

# ===============================
# PROJECT ROOT DIRECTORY
# ===============================
# This is the base directory where all inputs and outputs will be stored.
# Change this path to your own project location.

PROJECT_BASE_DIR="/path/to/your/project"

# ===============================
# INPUT DATA
# ===============================

# Directory containing paired-end FASTQ files (Phred33)
INPUT_FASTQ_DIR="${PROJECT_BASE_DIR}/raw_fastq"

# Sample metadata file (QIIME2-compatible TSV)
SAMPLE_METADATA_FILE="${PROJECT_BASE_DIR}/sample-metadata.tsv"

# Manifest file for importing FASTQ files into QIIME2
MANIFEST_FILE="${PROJECT_BASE_DIR}/manifest.tsv"

# ===============================
# TAXONOMIC CLASSIFIER
# ===============================
#
# IMPORTANT:
# The classifier MUST match:
#   - Amplicon region (e.g. V3–V4)
#   - Read length after trimming/truncation
#
# Using a generic classifier is a common source of poor taxonomic resolution.

SILVA_CLASSIFIER_QZA="/path/to/silva-V3V4-custom-classifier.qza"

# ===============================
# COMPUTATIONAL RESOURCES
# ===============================
#
# Adjust according to your system.
# High values improve speed but increase memory usage.

N_THREADS=16

###############################################################################
# OUTPUT DIRECTORIES
###############################################################################
#
# All outputs are organized into subdirectories for clarity.
# These directories will be created automatically if they do not exist.
#
###############################################################################

QIIME2_OUTPUT_DIR="${PROJECT_BASE_DIR}/qiime2_outputs"
VISUALIZATION_DIR="${PROJECT_BASE_DIR}/qiime2_visualizations"
EXPORT_DIR="${PROJECT_BASE_DIR}/exported_data"

mkdir -p "${QIIME2_OUTPUT_DIR}"
mkdir -p "${VISUALIZATION_DIR}"
mkdir -p "${EXPORT_DIR}"

###############################################################################
# END OF CONFIGURATION BLOCK
###############################################################################


###############################################################################
# STEP 1 — Import raw sequencing data into QIIME2
###############################################################################
#
# Raw FASTQ files are imported using a manifest file.
# This step does NOT modify the data; it only makes them readable by QIIME2.
#
# It is critical that:
# - Phred encoding is correct
# - Paired-end reads are correctly labeled
#
###############################################################################

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "${MANIFEST_FILE}" \
  --output-path "${QIIME2_OUTPUT_DIR}/demux-paired-end.qza" \
  --input-format PairedEndFastqManifestPhred33V2

###############################################################################
# STEP 2 — Visualize sequencing quality
###############################################################################
#
# This visualization is ESSENTIAL.
# All downstream parameter choices (truncation, trimming) depend on it.
#
###############################################################################

qiime demux summarize \
  --i-data "${QIIME2_OUTPUT_DIR}/demux-paired-end.qza" \
  --o-visualization "${VISUALIZATION_DIR}/demux-summary.qzv"

###############################################################################
# STEP 3 — Primer trimming with cutadapt
###############################################################################
#
# Even if primers were removed by the sequencing provider,
# trimming them again ensures:
# - Exact removal
# - Consistent read starts
# - Improved denoising performance
#
# IMPORTANT:
# Replace primer sequences with those used in your experiment.
#
###############################################################################

FORWARD_PRIMER="CCTACGGGNGGCWGCAG"
REVERSE_PRIMER="GACTACHVGGGTATCTAATCC"

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences "${QIIME2_OUTPUT_DIR}/demux-paired-end.qza" \
  --p-front-f "${FORWARD_PRIMER}" \
  --p-front-r "${REVERSE_PRIMER}" \
  --p-discard-untrimmed \
  --p-cores "${N_THREADS}" \
  --o-trimmed-sequences "${QIIME2_OUTPUT_DIR}/demux-trimmed.qza"

###############################################################################
# STEP 4 — Quality check after primer trimming
###############################################################################
#
# This step confirms that trimming behaved as expected.
#
###############################################################################

qiime demux summarize \
  --i-data "${QIIME2_OUTPUT_DIR}/demux-trimmed.qza" \
  --o-visualization "${VISUALIZATION_DIR}/demux-trimmed-summary.qzv"

###############################################################################
# END OF STEP 1
###############################################################################

###############################################################################
# STEP 5 — Denoising and ASV inference with DADA2
###############################################################################
#
# DADA2 performs:
# - Quality filtering
# - Error model learning
# - Dereplication
# - Chimera removal
#
# The result is a table of Amplicon Sequence Variants (ASVs),
# which represent exact biological sequences inferred from the data.
#
# IMPORTANT CONCEPT:
# ASVs are NOT OTUs.
# They are inferred sequences with single-nucleotide resolution.
#
###############################################################################
#
# PARAMETER SELECTION
#
# The truncation lengths below MUST be chosen after inspecting:
#   - demux-summary.qzv
#   - demux-trimmed-summary.qzv
#
# There is NO universally correct value.
#
###############################################################################

# Example truncation values (to be adjusted by the user)
TRUNC_LEN_F=240
TRUNC_LEN_R=200

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "${QIIME2_OUTPUT_DIR}/demux-trimmed.qza" \
  --p-trunc-len-f "${TRUNC_LEN_F}" \
  --p-trunc-len-r "${TRUNC_LEN_R}" \
  --p-n-threads "${N_THREADS}" \
  --o-table "${QIIME2_OUTPUT_DIR}/table-dada2.qza" \
  --o-representative-sequences "${QIIME2_OUTPUT_DIR}/rep-seqs-dada2.qza" \
  --o-denoising-stats "${QIIME2_OUTPUT_DIR}/denoising-stats.qza"

###############################################################################
# STEP 6 — DADA2 summary statistics
###############################################################################


###############################################################################
# STEP 7 — Taxonomic classification of ASVs
###############################################################################
#
# Taxonomic assignment is performed using a Naive Bayes classifier.
#
# ⚠️ CRITICAL CONCEPT:
# The classifier MUST be consistent with:
#   - The amplified region (e.g. V3–V4)
#   - The final ASV length after trimming/truncation
#
# Using a generic full-length classifier is a common source of:
#   - Low confidence assignments
#   - Poor genus-level resolution
#   - Misleading biological interpretations
#
###############################################################################
#
# CLASSIFIER STRATEGY USED IN THIS WORKFLOW
#
# - Custom-trained classifier
# - Region-specific (V3–V4)
# - Length-matched to ASVs
# - Curated using RESCRIPt to ensure:
#     - Reduced redundancy
#     - Lower memory usage
#     - Pipeline stability
#
# IMPORTANT:
# RESCRIPt does NOT artificially increase taxonomic resolution.
# Its primary value is methodological robustness and reproducibility.
#
###############################################################################

qiime feature-classifier classify-sklearn \
  --i-classifier "${SILVA_CLASSIFIER_QZA}" \
  --i-reads "${QIIME2_OUTPUT_DIR}/rep-seqs-dada2.qza" \
  --p-n-jobs "${N_THREADS}" \
  --o-classification "${QIIME2_OUTPUT_DIR}/taxonomy.qza"

###############################################################################
# STEP 8 — Taxonomy visualization
###############################################################################
#
# This table reports confidence scores and taxonomic depth.
# Users are strongly encouraged to inspect:
#   - Mean confidence
#   - Fraction of ASVs classified to genus level
#
###############################################################################

qiime metadata tabulate \
  --m-input-file "${QIIME2_OUTPUT_DIR}/taxonomy.qza" \
  --o-visualization "${VISUALIZATION_DIR}/taxonomy.qzv"

###############################################################################
# END OF TAXONOMIC CLASSIFICATION
###############################################################################

###############################################################################
# NOTE ON CLASSIFIER TRAINING AND RESCRIPt
###############################################################################
#
# Classifier training is intentionally NOT included in this script.
#
# Training a classifier is:
#   - Computationally expensive
#   - Dependent on reference database choice
#   - Dependent on available hardware
#
# In this workflow, classifier training and validation are:
#   - Documented separately
#   - Versioned
#   - Justified
#
# See:
#   02_qiime2_pipeline/classifiers/README_classifiers.md
#
# RESCRIPt was used to:
#   - Remove low-quality reference sequences
#   - Dereplicate redundant entries
#   - Enable classifier training on limited hardware
#
# Importantly:
#   - Taxonomic assignment gains were driven by length consistency
#   - Not by aggressive database pruning
#
###############################################################################
