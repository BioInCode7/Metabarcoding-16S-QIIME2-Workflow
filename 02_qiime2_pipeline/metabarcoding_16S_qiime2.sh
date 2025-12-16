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

# ===============================
# PIPELINE CONTROL FLAGS
# ===============================
# Set these to 1 to stop the pipeline at a given step

STOP_AFTER_MULTIQC=0
STOP_AFTER_DADA2=0

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
# Resolve repository root directory (script lives in 02_qiime2_pipeline/)
PROJECT_BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ===============================
# INPUT DATA
# ===============================

# Directory containing paired-end FASTQ files (Phred33)
INPUT_FASTQ_DIR="${PROJECT_BASE_DIR}/01_input_example/fastq_subsampled"

# Sample metadata file (QIIME2-compatible TSV)
SAMPLE_METADATA_FILE="${PROJECT_BASE_DIR}/01_input_example/sample-metadata.txt"

# Manifest file for importing FASTQ files into QIIME2
MANIFEST_FILE="${PROJECT_BASE_DIR}/01_input_example/manifest_v2.tsv"

# DEBUG (ahora sí)
echo "DEBUG: PROJECT_BASE_DIR = ${PROJECT_BASE_DIR}"
echo "DEBUG: INPUT_FASTQ_DIR = ${INPUT_FASTQ_DIR}"
ls -lh "${INPUT_FASTQ_DIR}"


###############################################################################
# INPUT VALIDATION
###############################################################################

echo "🔍 Validating input files..."

[[ -d "$INPUT_FASTQ_DIR" ]] || {
  echo "ERROR: FASTQ directory not found: $INPUT_FASTQ_DIR"
  exit 1
}

[[ -f "$MANIFEST_FILE" ]] || {
  echo "ERROR: Manifest file not found: $MANIFEST_FILE"
  exit 1
}

[[ -f "$SAMPLE_METADATA_FILE" ]] || {
  echo "ERROR: Sample metadata file not found: $SAMPLE_METADATA_FILE"
  exit 1
}

echo "✅ Input validation passed"

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

SILVA_CLASSIFIER_QZA="${PROJECT_BASE_DIR}/02_qiime2_pipeline/classifiers/silva-138-v3v4-420bp/silva-138-V3V4-420-classifier.qza"

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
# STEP 0A — Quality control of raw FASTQ files (FastQC)
###############################################################################
#
# FastQC provides per-base quality scores, GC content, adapter contamination,
# and overrepresented sequence diagnostics.
#
# This step is CRITICAL for:
#   - Deciding trimming and truncation parameters
#   - Assessing sequencing quality BEFORE any processing
#
# FastQC is run on the raw (subsampled) FASTQ files included in the repository.
#
###############################################################################

FASTQC_DIR="${PROJECT_BASE_DIR}/qc_fastqc"
MULTIQC_DIR="${PROJECT_BASE_DIR}/qc_multiqc"

mkdir -p "${FASTQC_DIR}"
mkdir -p "${MULTIQC_DIR}"

echo "▶ Running FastQC on raw FASTQ files..."

fastqc \
  --threads "${N_THREADS}" \
  --outdir "${FASTQC_DIR}" \
  "${INPUT_FASTQ_DIR}"/*.fastq.gz

echo "✅ FastQC completed"

###############################################################################
# STEP 0B — Aggregate QC reports with MultiQC
###############################################################################
#
# MultiQC aggregates all FastQC reports into a single interactive HTML file.
#
# This report is the PRIMARY document used to:
#   - Decide truncation lengths for DADA2
#   - Evaluate overall sequencing quality
#
# This file should be archived and cited in reports if relevant.
#
###############################################################################

echo "▶ Running MultiQC..."

multiqc \
  "${FASTQC_DIR}" \
  --outdir "${MULTIQC_DIR}"

echo "✅ MultiQC report generated:"
echo "   ${MULTIQC_DIR}/multiqc_report.html"

if [[ "${STOP_AFTER_MULTIQC}" -eq 1 ]]; then
  echo "🛑 Pipeline stopped after MultiQC as requested."
  exit 0
fi


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
TRUNC_LEN_F=280
TRUNC_LEN_R=260

# NOTE:
# Conservative truncation values were initially tested (270/240),
# but inspection of demux summary plots showed consistently high
# quality scores beyond these positions. Therefore, truncation
# lengths were extended to retain additional high-quality bases
# and maximize taxonomic resolution.

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "${QIIME2_OUTPUT_DIR}/demux-trimmed.qza" \
  --p-trunc-len-f "${TRUNC_LEN_F}" \
  --p-trunc-len-r "${TRUNC_LEN_R}" \
  --p-n-threads "${N_THREADS}" \
  --o-table "${QIIME2_OUTPUT_DIR}/table-dada2.qza" \
  --o-representative-sequences "${QIIME2_OUTPUT_DIR}/rep-seqs-dada2.qza" \
  --o-denoising-stats "${QIIME2_OUTPUT_DIR}/denoising-stats.qza" \
  --o-base-transition-stats "${QIIME2_OUTPUT_DIR}/base-transition-stats.qza"


  if [[ "${STOP_AFTER_DADA2}" -eq 1 ]]; then
  echo "🛑 Pipeline stopped after DADA2 as requested."
  exit 0
fi

###############################################################################
# STEP 6 — DADA2 summary statistics
###############################################################################


# This step generates human-readable summaries of the DADA2 outputs.
# These visualizations are ESSENTIAL to:
#   - Evaluate read retention
#   - Detect problematic samples
#   - Justify denoising parameters
#
###############################################################################

# --------------------------------
# 6.1 Feature table summary
# --------------------------------
#
# Shows:
# - Number of features (ASVs) per sample
# - Total frequency per sample
#
# Used to:
# - Detect low-depth samples
# - Inform rarefaction depth choice
#

qiime feature-table summarize \
  --i-table "${QIIME2_OUTPUT_DIR}/table-dada2.qza" \
  --o-visualization "${VISUALIZATION_DIR}/table-dada2-summary.qzv" \
  --m-sample-metadata-file "${SAMPLE_METADATA_FILE}"

# --------------------------------
# 6.2 ASV sequence length distribution
# --------------------------------
#
# Confirms:
# - Consistency of ASV lengths
# - Absence of truncated or anomalous sequences
#

qiime feature-table tabulate-seqs \
  --i-data "${QIIME2_OUTPUT_DIR}/rep-seqs-dada2.qza" \
  --o-visualization "${VISUALIZATION_DIR}/rep-seqs-dada2.qzv"

# --------------------------------
# 6.3 DADA2 denoising statistics
# --------------------------------
#
# Reports:
# - Reads input
# - Reads filtered
# - Reads merged
# - Non-chimeric reads
#
# This is the PRIMARY diagnostic for DADA2 performance.
#

qiime metadata tabulate \
  --m-input-file "${QIIME2_OUTPUT_DIR}/denoising-stats.qza" \
  --o-visualization "${VISUALIZATION_DIR}/denoising-stats.qzv"

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
# STEP 9 — Taxonomic filtering of non-bacterial sequences
###############################################################################

qiime taxa filter-table \
  --i-table "${QIIME2_OUTPUT_DIR}/table-dada2.qza" \
  --i-taxonomy "${QIIME2_OUTPUT_DIR}/taxonomy.qza" \
  --p-exclude mitochondria,chloroplast,eukaryota \
  --o-filtered-table "${QIIME2_OUTPUT_DIR}/table-taxa-filtered.qza"

qiime feature-table summarize \
  --i-table "${QIIME2_OUTPUT_DIR}/table-taxa-filtered.qza" \
  --o-visualization "${VISUALIZATION_DIR}/table-taxa-filtered-summary.qzv" \
  --m-sample-metadata-file "${SAMPLE_METADATA_FILE}"

  ###############################################################################
# STEP 10 — Rarefaction
###############################################################################

SAMPLING_DEPTH=10000

qiime feature-table rarefy \
  --i-table "${QIIME2_OUTPUT_DIR}/table-taxa-filtered.qza" \
  --p-sampling-depth "${SAMPLING_DEPTH}" \
  --o-rarefied-table "${QIIME2_OUTPUT_DIR}/table-rarefied.qza"

  ###############################################################################
# STEP 11 — Phylogenetic tree construction
###############################################################################

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${QIIME2_OUTPUT_DIR}/rep-seqs-dada2.qza" \
  --o-alignment "${QIIME2_OUTPUT_DIR}/aligned-rep-seqs.qza" \
  --o-masked-alignment "${QIIME2_OUTPUT_DIR}/masked-aligned-rep-seqs.qza" \
  --o-tree "${QIIME2_OUTPUT_DIR}/unrooted-tree.qza" \
  --o-rooted-tree "${QIIME2_OUTPUT_DIR}/rooted-tree.qza"

  ###############################################################################
# STEP 12 — Core diversity metrics
###############################################################################

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "${QIIME2_OUTPUT_DIR}/rooted-tree.qza" \
  --i-table "${QIIME2_OUTPUT_DIR}/table-rarefied.qza" \
  --p-sampling-depth "${SAMPLING_DEPTH}" \
  --m-metadata-file "${SAMPLE_METADATA_FILE}" \
  --output-dir "${QIIME2_OUTPUT_DIR}/core-metrics"

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
###############################################################################
# STEP V — Beta diversity statistical testing (PERMANOVA)
###############################################################################
#
# PERMANOVA (Permutational Multivariate Analysis of Variance) is used to test
# whether the centroids of predefined groups differ in multivariate space.
#
# It operates on distance matrices derived from beta diversity metrics.
#
###############################################################################
#
# WHAT PERMANOVA TESTS
#
# - Null hypothesis:
#     The centroids of the groups are identical in multivariate space.
#
# - Alternative hypothesis:
#     At least one group centroid differs.
#
###############################################################################
#
# WHAT PERMANOVA DOES NOT TEST
#
# - It does NOT identify which taxa drive differences
# - It does NOT measure within-group dispersion
# - It does NOT imply causality
#
# A significant PERMANOVA result indicates group separation,
# NOT mechanistic explanation.
#
###############################################################################
#
# IMPORTANT ASSUMPTION
#
# PERMANOVA assumes homogeneity of multivariate dispersion.
#
# Therefore:
# - Significant results SHOULD be accompanied by a dispersion test
# - Differences in dispersion can inflate false positives
#
###############################################################################
#
# CHOICE OF METADATA VARIABLE
#
# The grouping variable must:
#   - Reflect the experimental design
#   - Be defined a priori
#   - Avoid circular definitions
#
###############################################################################

###############################################################################
# PERMANOVA ON BRAY–CURTIS DISTANCE
###############################################################################

qiime diversity beta-group-significance \
  --i-distance-matrix "${QIIME2_OUTPUT_DIR}/core-metrics/bray_curtis_distance_matrix.qza" \
  --m-metadata-file "${SAMPLE_METADATA_FILE}" \
  --m-metadata-column group \
  --p-method permanova \
  --p-permutations 999 \
  --o-visualization "${VISUALIZATION_DIR}/permanova-bray-curtis-group.qzv"

###############################################################################
# PERMANOVA ON WEIGHTED UNIFRAC DISTANCE
###############################################################################
#
# Weighted UniFrac incorporates both:
#   - Relative abundance
#   - Phylogenetic relationships
#
###############################################################################

qiime diversity beta-group-significance \
  --i-distance-matrix "${QIIME2_OUTPUT_DIR}/core-metrics/weighted_unifrac_distance_matrix.qza" \
  --m-metadata-file "${SAMPLE_METADATA_FILE}" \
  --m-metadata-column group \
  --p-method permanova \
  --p-permutations 999 \
  --o-visualization "${VISUALIZATION_DIR}/permanova-weighted-unifrac-group.qzv"

###############################################################################
# DISPERSION TEST (BETADISPER)
###############################################################################
#
# This test evaluates whether group dispersions differ.
# A significant dispersion test requires cautious interpretation
# of PERMANOVA results.
#
###############################################################################

qiime diversity beta-group-significance \
  --i-distance-matrix "${QIIME2_OUTPUT_DIR}/core-metrics/bray_curtis_distance_matrix.qza" \
  --m-metadata-file "${SAMPLE_METADATA_FILE}" \
  --m-metadata-column group \
  --p-method permdisp \
  --p-permutations 999 \
  --o-visualization "${VISUALIZATION_DIR}/permdisp-bray-curtis-group.qzv"

###############################################################################
# END OF PERMANOVA ANALYSES
###############################################################################

###############################################################################
# NOTE ON REPORTING RESULTS (PAPER-READY)
###############################################################################
#
# Example phrasing:
#
# "Community composition differed between groups based on Bray–Curtis
# dissimilarities (PERMANOVA, 999 permutations, p < 0.05). However, dispersion
# tests indicated [no / significant] differences in within-group variability,
# and results were interpreted accordingly."
#
###############################################################################


###############################################################################
# STEP U — Functional inference with PICRUSt2
###############################################################################
#
# PICRUSt2 predicts the FUNCTIONAL POTENTIAL of microbial communities
# based on 16S rRNA gene sequences and reference genomes.
#
###############################################################################
#
# ⚠️ CRITICAL CONCEPTUAL WARNING
#
# PICRUSt2:
#   - Does NOT measure gene presence directly
#   - Does NOT measure gene expression
#   - Does NOT measure enzymatic activity
#
# It infers potential functions based on phylogenetic proximity
# to reference genomes.
#
# Therefore:
#   Predicted function ≠ realized function
#
###############################################################################
#
# WHAT PICRUSt2 CAN BE USED FOR
#
# - Hypothesis generation
# - Functional comparison between groups
# - Identification of overrepresented pathways
#
###############################################################################
#
# WHAT PICRUSt2 CANNOT PROVE
#
# - Active metabolism
# - Biodegradation rates
# - Enzymatic efficiency
# - Causality between taxa and phenotype
#
###############################################################################
#
# RECOMMENDED INTERPRETATION
#
# PICRUSt2 results should be interpreted as:
#   - Functional potential
#   - Community-level tendencies
#   - Supporting evidence, never as standalone proof
#
###############################################################################
#
# IMPORTANT METHODOLOGICAL NOTE
#
# PICRUSt2 MUST be run on:
#   - Non-rarefied feature table
#   - ASV representative sequences
#
# Rarefaction removes quantitative information required
# for functional inference.
#
###############################################################################


###############################################################################
# EXPORT DATA FOR PICRUSt2
###############################################################################
#
# PICRUSt2 is run outside QIIME2.
# Therefore, data must be exported to standard formats.
#
###############################################################################

# Export feature table (BIOM format)
qiime tools export \
  --input-path "${QIIME2_OUTPUT_DIR}/table-taxa-filtered.qza" \
  --output-path "${EXPORT_DIR}/picrust2"

# Export representative sequences (FASTA)
qiime tools export \
  --input-path "${QIIME2_OUTPUT_DIR}/rep-seqs-dada2.qza" \
  --output-path "${EXPORT_DIR}/picrust2"

###############################################################################
# Convert BIOM to TSV (optional but recommended for inspection)
###############################################################################

biom convert \
  -i "${EXPORT_DIR}/picrust2/feature-table.biom" \
  -o "${EXPORT_DIR}/picrust2/feature-table.tsv" \
  --to-tsv

###############################################################################
# END OF EXPORTS FOR PICRUSt2
###############################################################################


###############################################################################
# NOTE ON RUNNING PICRUSt2
###############################################################################
#
# PICRUSt2 should be executed in a dedicated environment.
#
# Example command (to be run outside this script):
#
# picrust2_pipeline.py \
#   -s rep-seqs-dada2.fasta \
#   -i feature-table.biom \
#   -o picrust2_out \
#   -p 8
#
# Downstream analyses (e.g. pathway comparison, visualization)
# are performed in R.
#
###############################################################################
#
# FINAL INTERPRETATION GUIDELINE (PAPER-READY)
#
# "Functional profiles were inferred using PICRUSt2 and interpreted as
# community-level functional potential. These predictions do not constitute
# direct evidence of gene presence, expression, or metabolic activity."
#
###############################################################################
#
# END OF PICRUSt2 SECTION
###############################################################################

