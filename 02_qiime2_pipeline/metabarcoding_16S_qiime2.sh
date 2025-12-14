#!/bin/bash
################################################################################
# METABARCODING 16S rRNA PIPELINE (QIIME2)
#
# PURPOSE
# -------
# This script implements a complete, reproducible, and documented workflow
# for 16S rRNA metabarcoding analysis using QIIME2 and downstream tools.
#
# The script is intentionally verbose and commented. It is meant to be read
# as a methodological document, not just executed.
#
# The example system (insect gut microbiota under plastic-based diet) is used
# ONLY as a case study. The pipeline is transferable to other systems.
#
# WHAT THIS PIPELINE DOES
# ----------------------
# - Processes raw paired-end FASTQ files
# - Performs quality control and denoising (DADA2)
# - Generates ASVs (Amplicon Sequence Variants)
# - Assigns taxonomy using a trained classifier
# - Computes alpha and beta diversity
# - Performs multivariate statistical testing (PERMANOVA)
# - Exports data for downstream analysis in R (e.g. PICRUSt, plots)
#
# WHAT THIS PIPELINE DOES NOT DO
# -----------------------------
# - It does NOT demonstrate biodegradation
# - It does NOT measure functional activity
# - It does NOT infer causality
#
# Functional profiles derived from 16S (e.g. PICRUSt) are INFERENCES ONLY.
#
# REPRODUCIBILITY
# ---------------
# The workflow is reproducible at the pipeline level.
# Results depend on sequencing depth, primers, classifiers, and design.
################################################################################
################################################################################
# 1. USER CONFIGURATION
#
# IMPORTANT:
# ----------
# Modify ONLY this section.
# The rest of the script assumes these variables are correct.
################################################################################

# Base directory of the project
PROJECT_DIR="/path/to/your/project"

# Directory containing raw FASTQ files
INPUT_FASTQ_DIR="/path/to/fastq"

# Sample metadata file (QIIME2-compliant TSV)
METADATA_FILE="/path/to/sample-metadata.tsv"

# Taxonomic classifier (.qza)
# CRITICAL NOTE:
# --------------
# The classifier MUST:
# - match your primer set
# - match your amplicon region
# - match your final ASV length
#
# Using a generic classifier is strongly discouraged.
CLASSIFIER_QZA="/path/to/classifier.qza"

# Primer sequences (leave empty if not trimming primers)
FORWARD_PRIMER=""
REVERSE_PRIMER=""

# Number of CPU threads
THREADS=16
################################################################################
# 2. DIRECTORY STRUCTURE
#
# A clear directory structure improves:
# - reproducibility
# - debugging
# - long-term project maintenance
################################################################################

mkdir -p "${PROJECT_DIR}"/{manifests,qiime_imported,qiime_cutadapt,qiime_dada2}
mkdir -p "${PROJECT_DIR}"/{qiime_downstream,qiime_visualizations}
mkdir -p "${PROJECT_DIR}"/{exported_data,summary_stats}
################################################################################
# 3. INPUT DATA AND MANIFEST
#
# QIIME2 does not read FASTQ files directly.
# Instead, it requires a MANIFEST file linking sample IDs to file paths.
#
# This design:
# - enforces explicit sample tracking
# - prevents silent mismatches
################################################################################

# The user must provide a manifest.tsv
# Example manifest is provided in the repository.
MANIFEST_FILE="${PROJECT_DIR}/manifests/manifest.tsv"
################################################################################
# 4. IMPORT RAW READS INTO QIIME2
#
# This step converts raw FASTQ files into a QIIME2 artifact (.qza).
# No processing occurs here.
################################################################################

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$MANIFEST_FILE" \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path "${PROJECT_DIR}/qiime_imported/demux_raw.qza"
################################################################################
# 5. ADAPTER AND PRIMER TRIMMING (CUTADAPT)
#
# WHY?
# ----
# Primer and adapter sequences:
# - are not biological information
# - bias error models
# - interfere with taxonomy assignment
#
# Removing them improves:
# - DADA2 denoising
# - taxonomic resolution
################################################################################

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences "${PROJECT_DIR}/qiime_imported/demux_raw.qza" \
  --p-front-f "$FORWARD_PRIMER" \
  --p-front-r "$REVERSE_PRIMER" \
  --p-cores "$THREADS" \
  --o-trimmed-sequences "${PROJECT_DIR}/qiime_cutadapt/demux_trimmed.qza"
################################################################################
# 6. DENOISING WITH DADA2
#
# DADA2 performs:
# - error modeling
# - dereplication
# - chimera removal
#
# Output:
# - ASVs (exact sequence variants, not OTUs)
#
# PARAMETER CHOICE IS CRITICAL.
################################################################################

# Example parameters (MUST be optimized per dataset)
TRUNC_LEN=400
TRIM_LEFT=0

qiime dada2 denoise-single \
  --i-demultiplexed-seqs "${PROJECT_DIR}/qiime_cutadapt/demux_trimmed.qza" \
  --p-trunc-len "$TRUNC_LEN" \
  --p-trim-left "$TRIM_LEFT" \
  --p-n-threads "$THREADS" \
  --o-table "${PROJECT_DIR}/qiime_dada2/table.qza" \
  --o-representative-sequences "${PROJECT_DIR}/qiime_dada2/rep_seqs.qza" \
  --o-denoising-stats "${PROJECT_DIR}/qiime_dada2/stats.qza"
################################################################################
# 7. TAXONOMIC CLASSIFICATION
#
# Taxonomy is assigned using a Naive Bayes classifier.
#
# IMPORTANT:
# ----------
# Classifier quality has a stronger effect on results
# than almost any other downstream step.
################################################################################

qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER_QZA" \
  --i-reads "${PROJECT_DIR}/qiime_dada2/rep_seqs.qza" \
  --p-n-jobs "$THREADS" \
  --o-classification "${PROJECT_DIR}/qiime_downstream/taxonomy.qza"
################################################################################
# 8. FILTERING AND DIVERSITY ANALYSES
#
# Includes:
# - contaminant removal
# - rarefaction
# - alpha diversity
# - beta diversity
#
# Rarefaction is used ONLY for diversity metrics.
################################################################################
################################################################################
# 9. EXPORT FOR DOWNSTREAM ANALYSIS
#
# Exports:
# - feature table
# - taxonomy
# - representative sequences
#
# These files can be used in:
# - R (ggplot, phyloseq)
# - PICRUSt
################################################################################
################################################################################
# FINAL NOTES
#
# - 16S provides relative abundances, not absolute counts.
# - Functional inference ≠ functional measurement.
# - Results must be interpreted in ecological context.
#
# This pipeline is a tool.
# Biological interpretation is the responsibility of the researcher.
################################################################################
