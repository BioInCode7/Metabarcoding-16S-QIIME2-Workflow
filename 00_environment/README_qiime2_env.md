# QIIME2 environment

This workflow was developed, tested and validated using:

- **QIIME2 Amplicon distribution 2025.10**
- Python 3.10
- Conda-based installation

## Recommended installation (official)

Users are strongly encouraged to install **the same QIIME2 release**
to ensure full compatibility with this workflow.

Follow the official QIIME2 installation instructions for the
Amplicon distribution:

https://docs.qiime2.org/2025.10/install/

For Linux (example):

```bash
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2025.10-py310-linux-conda.yml
conda env create -n qiime2-amplicon-2025.10 --file qiime2-amplicon-2025.10-py310-linux-conda.yml
conda activate qiime2-amplicon-2025.10

