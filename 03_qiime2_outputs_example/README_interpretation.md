# QIIME2 results – Quality control and denoising

This document summarizes the interpretation of the main QIIME2 visualizations
generated during the 16S metabarcoding pipeline.

⚠️ Note: The analyses were performed on subsampled datasets (20,000 reads per sample)
for workflow testing and reproducibility. Results are not intended for biological inference.

## 1. Demultiplexing summary (`demux-summary.qzv`)

- All samples show successful import of paired-end reads.
- Forward and reverse reads have a consistent length of ~301 bp.
- Quality profiles indicate:
  - High median quality scores (>Q30) across most of the read length.
  - A gradual quality decrease towards the 3' end, more pronounced in reverse reads.

Conclusion:
The sequencing quality is sufficient for denoising with DADA2, with appropriate
truncation parameters applied in the next step.

## 2. Trimmed reads summary (`demux-trimmed-summary.qzv`)

- After primer trimming, read lengths are reduced to:
  - Forward reads: median ~282 bp
  - Reverse reads: median ~278 bp
- Quality profiles remain stable after trimming.
- No abnormal quality drops or length heterogeneity are observed.

Conclusion:
Primer removal was successful and did not introduce artifacts or excessive read loss.
The trimmed reads are suitable for DADA2 denoising and merging.

## 3. DADA2 denoising

DADA2 successfully generated:

- A feature table (`table-dada2.qza`)
- Representative sequences (`rep-seqs-dada2.qza`)
- Denoising statistics (`denoising-stats.qza`)

No errors were observed during denoising, indicating that trimming and truncation
parameters were appropriate for the quality profiles observed.

Detailed statistics will be summarized after taxonomy assignment and filtering.

