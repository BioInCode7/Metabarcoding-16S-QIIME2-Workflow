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
