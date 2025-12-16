# Metabarcoding-16S-QIIME2-Workflow

A **reproducible, conservative and well-documented workflow** for 16S rRNA metabarcoding analysis using **QIIME2**, with downstream analysis in **R** and functional inference using **PICRUSt2**.

This repository is designed as a **methodological guide**, not as a data repository.

---

## 🎯 Scope and philosophy

This workflow aims to:

- Provide a **single, fully commented QIIME2 script** covering the complete 16S metabarcoding pipeline
- Explain **why each tool and parameter is used**
- Promote **honest interpretation** of metabarcoding results
- Serve as a **template adaptable to any biological system**

> ⚠️ This repository does **not** claim experimental validation of microbial functions.  
> Functional profiles inferred from 16S data (e.g. via PICRUSt2) represent **predicted potential**, not measured activity.

---

## ❌ What this repository does NOT do

This workflow **does not**:

- Demonstrate biodegradation, metabolic activity, or functional validation
- Replace shotgun metagenomics or experimental assays
- Provide biological conclusions beyond what 16S data can support

> **Inference is not evidence of function.**

---

## 🧬 Example system (illustrative only)

The example structure and scripts were developed using a real case study involving:

- Insect gut microbiota
- Extreme dietary conditions (e.g. plastic-based diet)
- Probiotic inoculation

⚠️ **Important**:  
The biological system is used **only as an example**.  
No raw sequencing data or real experimental outputs are provided.

Users are expected to:
- Supply their own FASTQ files
- Adapt metadata and experimental design to their system

For reproducibility and computational accessibility, example analyses in this
repository are performed on uniformly subsampled FASTQ files. This does not
replace full-depth analyses of the original sequencing data.


---

## 🧱 Repository structure

Metabarcoding-16S-QIIME2-Workflow/
├── 00_environment/ # Conda environments (QIIME2, R)
├── 01_input_example/ # Example metadata and manifests (no real data)
├── 02_qiime2_pipeline/ # Core QIIME2 workflow (single main script)
│ ├── metabarcoding_16S_qiime2.sh
│ └── classifiers/
├── 03_qiime2_outputs_example/ # Explanation of expected QIIME2 outputs
├── 04_downstream_R/ # R scripts (one script per figure)
├── 05_picrust_pipeline/ # PICRUSt2 workflow
├── 06_methodological_notes/ # Limits, pitfalls and interpretation notes
└── figures_example/ # Example figures (illustrative only)


---

## 🛠️ Tools used

- **QIIME2** (denoising, taxonomy, diversity)
- **DADA2**
- **SILVA 138** reference database
- **RESCRIPt** (classifier curation)
- **PICRUSt2** (functional inference)
- **R** (ggplot2, phyloseq, vegan, tidyverse)

---

## 🔁 Reproducibility

This repository prioritizes:

- Reproducible **pipelines**, not reproducible **results**
- Explicit documentation of parameters
- Clear separation between raw data, processing and interpretation

---

## 📌 Intended audience

- Researchers new to metabarcoding
- Biologists transitioning into bioinformatics
- Anyone seeking a **transparent, non-inflated** 16S workflow

---

## 📖 Citation

If you use or adapt this workflow in academic work, please cite this repository.  
A `CITATION.cff` file is provided.

---

## 📬 Contact

Maintained by **Jesús Salinas**  
University of Almería  
LinkedIn: https://www.linkedin.com/in/jsalinasbiotech/