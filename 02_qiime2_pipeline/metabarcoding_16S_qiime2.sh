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

PROJECT_BASE_DIR="${PROJECT_BASE_DIR:-$(pwd)}"

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
# STEP 0 — Auto-generate paired-end FASTQ manifest (example data)
###############################################################################
#
# This pipeline automatically generates a QIIME2-compatible manifest file
# from the subsampled FASTQ files included in the repository.
#
# Expected filename pattern:
#   <sample-id>_R1_*.fastq.gz
#   <sample-id>_R2_*.fastq.gz
#
###############################################################################

INPUT_FASTQ_DIR="${PROJECT_BASE_DIR}/raw_fastq"
MANIFEST_FILE="${PROJECT_BASE_DIR}/manifest.tsv"

mkdir -p "$(dirname "$MANIFEST_FILE")"

echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" \
  > "$MANIFEST_FILE"

for R1 in "${INPUT_FASTQ_DIR}"/*_R1_*.fastq.gz; do
  SAMPLE_ID=$(basename "$R1" | sed 's/_R1_.*\.fastq\.gz//')

  R2="${INPUT_FASTQ_DIR}/${SAMPLE_ID}_R2_subsampled.fastq.gz"

  if [[ ! -f "$R2" ]]; then
    echo "ERROR: Missing R2 file for sample ${SAMPLE_ID}"
    exit 1
  fi

  echo -e "${SAMPLE_ID}\t${R1}\t${R2}" \
    >> "$MANIFEST_FILE"
done

echo "--- Paired-end FASTQ manifest generated ---"
echo "    $MANIFEST_FILE"

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

###############################################################################
# STEP X — Rarefaction and Sampling Depth Selection
###############################################################################
#
# Rarefaction is applied exclusively for diversity analyses (alpha and beta).
# It is NOT used for:
#   - Taxonomic composition
#   - Differential abundance analysis
#   - Functional inference (PICRUSt2)
#
# Rarefaction enforces an equal sequencing depth across samples to allow
# meaningful diversity comparisons.
#
# ⚠️ IMPORTANT:
# Rarefaction ALWAYS discards data.
# The goal is NOT to avoid data loss, but to make it explicit and justified.
#
###############################################################################
#
# CONCEPTUAL FRAMEWORK
#
# Sampling depth is the number of sequences per sample retained after
# subsampling. It represents a trade-off between:
#
#   - Retaining all biological replicates
#   - Preserving as much within-sample diversity as possible
#
# There is no universally "correct" sampling depth.
# The chosen value must be justified based on the dataset.
#
###############################################################################
#
# HOW TO SELECT SAMPLING DEPTH (REQUIRED PROCEDURE)
#
# 1. Inspect the feature table AFTER:
#    - DADA2 denoising
#    - Taxonomic filtering (mitochondria, chloroplasts, eukaryotes removed)
#
#    Command:
#
#    qiime feature-table summarize \
#      --i-table table-taxa-filtered.qza \
#      --o-visualization table-summary-taxa-filtered.qzv \
#      --m-sample-metadata-file sample-metadata.tsv
#
# 2. Open the visualization:
#      qiime view table-summary-taxa-filtered.qzv
#
# 3. Identify:
#    - The minimum sequencing depth
#    - Samples or groups with systematically lower depth
#
###############################################################################
#
# RAREFACTION CURVES AS DECISION SUPPORT
#
# Rarefaction curves MUST be interpreted metric by metric:
#
#   - Shannon diversity:
#       * Typically reaches saturation at relatively low depths
#       * Robust to subsampling
#       * Recommended for interpretation when sequencing depth is limited
#
#   - Observed Features:
#       * Sensitive to sequencing depth
#       * May or may not reach saturation
#
#   - Chao1:
#       * Estimates unseen richness
#       * Often does NOT plateau
#       * Lack of saturation does NOT invalidate the analysis
#
# Failure of Chao1 to plateau indicates incomplete richness capture,
# NOT methodological error.
#
###############################################################################
#
# FINAL DECISION CRITERIA
#
# The selected sampling depth should:
#   - Retain all biological replicates whenever possible
#   - Fall below the minimum sample depth
#   - Be supported by Shannon rarefaction curves approaching saturation
#
# Richness-based metrics should be interpreted conservatively when curves
# do not plateau.
#
###############################################################################
#
# EXAMPLE STATEMENT (PAPER-READY)
#
# "The sampling depth was selected as a compromise between retaining all samples
# and preserving sequencing depth. Rarefaction curves indicated that Shannon
# diversity approached saturation at this depth, while richness estimators did
# not fully plateau, suggesting underestimation of absolute richness but
# allowing robust comparative diversity analyses."
#
###############################################################################
#
# RAREFY FEATURE TABLE
#
###############################################################################

qiime feature-table rarefy \
  --i-table table-taxa-filtered.qza \
  --p-sampling-depth 3000 \
  --o-rarefied-table table-rarefied.qza

###############################################################################
# END OF RAREFACTION STEP
###############################################################################


###############################################################################
# STEP Y — Phylogenetic tree construction
###############################################################################
#
# A phylogenetic tree is required for phylogeny-based diversity metrics
# such as Faith's PD and UniFrac distances.
#
# IMPORTANT:
# The tree is built from representative ASV sequences inferred by DADA2.
#
# This tree:
#   - Represents relationships between ASVs
#   - Does NOT represent an organismal phylogeny
#
###############################################################################
#
# TREE CONSTRUCTION PIPELINE
#
# 1. Multiple sequence alignment (MAFFT)
# 2. Masking of hypervariable positions
# 3. Tree inference (FastTree)
# 4. Rooting
#
###############################################################################

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${QIIME2_OUTPUT_DIR}/rep-seqs-dada2.qza" \
  --o-alignment "${QIIME2_OUTPUT_DIR}/aligned-rep-seqs.qza" \
  --o-masked-alignment "${QIIME2_OUTPUT_DIR}/masked-aligned-rep-seqs.qza" \
  --o-tree "${QIIME2_OUTPUT_DIR}/unrooted-tree.qza" \
  --o-rooted-tree "${QIIME2_OUTPUT_DIR}/rooted-tree.qza"

###############################################################################
# STEP Z — Alpha and Beta diversity analyses
###############################################################################
#
# Diversity analyses are performed on the rarefied feature table ONLY.
#
# This ensures:
#   - Equal sampling depth across samples
#   - Comparable diversity estimates
#
###############################################################################
#
# DIVERSITY METRICS USED
#
# Alpha diversity:
#   - Shannon: accounts for richness and evenness
#   - Observed Features: raw ASV richness
#   - Faith's PD: phylogenetic diversity
#
# Beta d


###############################################################################
# STEP W — Taxonomic filtering of non-bacterial sequences
###############################################################################
#
# 16S rRNA primers frequently amplify:
#   - Mitochondrial sequences
#   - Chloroplast sequences
#   - Occasionally eukaryotic contaminants
#
# These sequences:
#   - Are biologically real
#   - BUT are not informative for bacterial community analyses
#
# Therefore, they must be removed BEFORE diversity analyses.
#
###############################################################################
#
# IMPORTANT:
# Taxonomic filtering affects:
#   - Feature table
#   - Downstream diversity metrics
#
# It does NOT affect:
#   - Raw ASV inference
#   - Reproducibility of the denoising step
#
###############################################################################

qiime taxa filter-table \
  --i-table "${QIIME2_OUTPUT_DIR}/table-dada2.qza" \
  --i-taxonomy "${QIIME2_OUTPUT_DIR}/taxonomy.qza" \
  --p-exclude mitochondria,chloroplast,eukaryota \
  --o-filtered-table "${QIIME2_OUTPUT_DIR}/table-taxa-filtered.qza"

###############################################################################
# STEP W2 — Summary after taxonomic filtering
###############################################################################
#
# This summary allows inspection of:
#   - Retained sequencing depth
#   - Potential sample loss
#
###############################################################################

qiime feature-table summarize \
  --i-table "${QIIME2_OUTPUT_DIR}/table-taxa-filtered.qza" \
  --o-visualization "${VISUALIZATION_DIR}/table-taxa-filtered-summary.qzv" \
  --m-sample-metadata-file "${SAMPLE_METADATA_FILE}"

###############################################################################
# END OF TAXONOMIC FILTERING
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

