# Taxonomic Classifier Training

This directory contains documentation and example scripts for training
custom Naive Bayes classifiers for 16S rRNA gene amplicon data.

Classifiers are trained to match:
- Amplicon region (V3–V4)
- Primer sequences
- ASV length distribution inferred by DADA2

---

## Reference Database

- Database: SILVA 138
- Taxonomy: 99% identity
- Source files:
  - silva-138-99-seqs.qza
  - silva-138-99-tax.qza

---

## Primer Sequences

The following primers were used for both sequencing and classifier training:

- Forward: CCTACGGGNGGCWGCAG
- Reverse: GACTACHVGGGTATCTAATCC

---

## Length Matching Strategy

Based on DADA2 results:
- Mean ASV length: ~407 bp
- Interquartile range: ~403–427 bp

Therefore, classifiers were trained using reference sequences truncated to
lengths close to the observed ASV distribution.

---

## Trained Classifiers

Classifiers were trained for multiple truncation lengths to evaluate performance:

- 400 bp
- 450 bp

Training was performed using the following steps:
1. Primer-based extraction of reference reads
2. Length truncation
3. Naive Bayes classifier training

The recommended classifier for this dataset is the 400 bp V3–V4 classifier,
which best matches the empirical ASV length distribution.

---

## Reproducibility

All classifier training steps are documented in:

train_classifier_example.sh


Users are encouraged to retrain classifiers if:
- Different primers are used
- Different truncation parameters are selected
- A different reference database is preferred
