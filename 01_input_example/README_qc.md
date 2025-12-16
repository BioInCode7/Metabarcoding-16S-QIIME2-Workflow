# Input example data and quality control (16S rRNA, V3–V4)

This directory contains **subsampled paired-end FASTQ files** together with
all quality control (QC) outputs used to demonstrate, validate, and document
the QIIME2 metabarcoding workflow implemented in this repository.

The focus of this directory is **methodological transparency**, not biological
interpretation.

---

## Important note on data usage

⚠️ **These FASTQ files do NOT correspond to full sequencing datasets.**

To ensure that the pipeline:

- runs quickly on standard hardware,
- is fully reproducible by third parties,
- and allows meaningful diversity comparisons,

all samples were **randomly subsampled to an equal sequencing depth**
prior to QIIME2 processing.

These data are intended exclusively for:
- pipeline development and testing,
- educational and methodological demonstration,
- example analyses accompanying this repository.

They are **not intended** to replace full-depth analyses of the original data.

---

## Subsampling strategy

- Tool: `seqtk`
- Reads retained per sample: **20,000 paired-end reads**
- Random seed: **42** (full reproducibility)
- Subsampling performed **before QIIME2 import**

This guarantees:
- identical sequencing depth across samples,
- valid alpha diversity comparisons,
- interpretable beta diversity and PERMANOVA results.

The subsampling procedure is fully documented in:
01_input_example/subsample_fastq.sh


Users analyzing their own data should either:
1. Skip the subsampling step, or
2. Adjust the subsampling depth according to their study design.

---

## FASTQ quality control (FastQC & MultiQC)

Raw subsampled FASTQ files were evaluated using **FastQC**, and reports were
aggregated with **MultiQC**.

### Summary of raw read quality

- Sequencing depth per sample: ~16,000–20,000 paired-end reads
- No evidence of adapter contamination
- No abnormal GC-content or duplication patterns
- Consistent quality profiles across all samples

MultiQC report:

01_input_example/qc_multiqc/multiqc_report.html


These results indicate high-quality Illumina sequencing data suitable for
conservative trimming and denoising.

---

## QIIME2 demultiplexing summary (`demux-summary.qzv`)

Raw reads were imported into QIIME2 using a paired-end manifest and summarized
with `qiime demux summarize`.

### Read length

- Forward reads: fixed length of **301 bp**
- Reverse reads: fixed length of **301 bp**

### Quality profiles

- Median quality scores > Q30 across most of the read length
- No abrupt quality drop-off prior to read ends
- High consistency across samples

These observations support retaining long read fragments while removing only
low-quality terminal regions.

---

## Primer trimming (`demux-trimmed-summary.qzv`)

Primers were removed using `qiime cutadapt trim-paired` to ensure:

- exact primer removal,
- consistent read starts,
- optimal DADA2 error modeling.

### Post-trimming read length distributions

**Forward reads**
- Median length: ~282 bp
- IQR: ~281–283 bp

**Reverse reads**
- Median length: ~278 bp
- IQR: ~277–279 bp

### Quality considerations

- Forward reads retained high quality (>Q30) across most positions
- Reverse reads showed a modest quality decline beyond ~280 bp
- Median reverse quality at position ~284 ≈ Q26

This pattern is typical of paired-end Illumina sequencing and informed
conservative truncation of reverse reads.

---

## DADA2 truncation parameter selection

Based on visual inspection of:

- `demux-summary.qzv`
- `demux-trimmed-summary.qzv`

the following truncation parameters were selected:

TRUNC_LEN_F = 280
TRUNC_LEN_R = 260


### Rationale

- Retains the majority of high-quality bases
- Removes low-quality tail regions (primarily in reverse reads)
- Ensures sufficient overlap for reliable paired-end merging
- Maximizes retained biological signal while maintaining denoising robustness

These parameters represent a data-driven compromise between read quality,
sequence length retention, and error-model stability.

---

## DADA2 denoising performance (`denoising-stats.qzv`)

DADA2 denoising showed consistent and stable performance across all samples:

- Reads passing quality filtering: ~85–87%
- Reads successfully merged: ~79–87%
- Non-chimeric reads retained: ~71–84%

No samples exhibited abnormal read loss or signs of poor error-model fitting.

---

## ASV characteristics

- Number of samples: 11
- Number of ASVs: ~370
- Mean ASV length: ~407 bp
- Length distribution consistent with expected V3–V4 amplicon size

These results confirm that the selected truncation parameters preserved
sufficient overlap and did not artificially shorten reconstructed sequences.

---

## Implications for taxonomic classifier training

Because DADA2 inference produced ASVs with a mean length of approximately
**407 bp**, taxonomic classifiers were trained using reference sequences
trimmed to comparable lengths.

This length matching avoids systematic biases associated with mismatched
classifier training and improves taxonomic assignment reliability.

Classifier training procedures are documented in:

02_qiime2_pipeline/classifiers/



---

## Manifest format

FASTQ files were imported using a **wide paired-end manifest** with one row
per sample and separate columns for forward and reverse reads:

- `sample-id`
- `forward-absolute-filepath`
- `reverse-absolute-filepath`

This format avoids duplicated sample IDs and is fully compatible with:

PairedEndFastqManifestPhred33V2


The manifest was generated programmatically using:

01_input_example/build_manifest.py


---

## Summary

All quality control steps and parameter choices in this directory are:

- empirically justified,
- fully documented,
- and reproducible.

This ensures transparency and methodological consistency throughout the
entire QIIME2 metabarcoding workflow.










