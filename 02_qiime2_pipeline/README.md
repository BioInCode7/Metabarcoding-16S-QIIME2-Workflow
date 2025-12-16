# QIIME2 Metabarcoding 16S Pipeline

This directory contains the main QIIME2-based metabarcoding workflow
for paired-end 16S rRNA amplicon data.

---

## RUNNING THE PIPELINE

This pipeline is executed as a Bash script and relies on a pre-existing
QIIME2 environment.

### Option 1 (recommended): define the project directory explicitly

Define the base directory of your project **before** running the script.
All input files and outputs will be resolved relative to this directory.

```bash
export PROJECT_BASE_DIR=/path/to/your/project
bash metabarcoding_16S_qiime2.sh

This option is recommended for reproducibility and clarity, especially
when running the pipeline from outside the project directory.

Option 2: run from within the project directory

If you run the script from inside your project directory, the pipeline
will automatically use the current working directory as the project base.

cd /path/to/your/project
bash metabarcoding_16S_qiime2.sh

Notes

The pipeline will create all required output directories automatically.

Input files (FASTQ, manifest, metadata) must already exist in the
project directory.

The value of PROJECT_BASE_DIR is resolved internally as:

PROJECT_BASE_DIR="${PROJECT_BASE_DIR:-$(pwd)}"
which allows both execution modes described above.

## Consistency Between DADA2 and Taxonomic Classification

DADA2 truncation parameters were selected empirically based on quality profile
inspection, resulting in ASVs with a mean length of approximately 407 bp.

To ensure methodological consistency, taxonomic classifiers were trained using
reference sequences trimmed to lengths matching the observed ASV distribution
(V3–V4 region, ~400–420 bp).

This avoids systematic biases associated with length-mismatched classifiers and
improves taxonomic assignment reliability.



