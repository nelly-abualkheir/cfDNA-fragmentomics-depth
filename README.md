# cfDNA Fragmentomics Depth Robustness Project

## Overview

This project evaluates how simulated reductions in sequencing depth affect the stability and diagnostic performance of cfDNA fragmentomic features for breast cancer classification.

The analysis compares global fragment-size features with genome-wide regional fragmentation patterns and evaluates classification performance as fragment counts are reduced from full depth to 20M, 10M, 5M, 2M, and 1M fragments.

## Research Question

How does simulated reduction in sequencing depth affect the stability of global and genome-wide cfDNA fragmentomic features, and which features maintain reliable cancer-versus-healthy classification at low fragment counts?

## Dataset

Public cfDNA whole-genome sequencing fragment data were obtained through FinaleDB from the Cristiano et al. (2019) study.

Final cohort:

- 46 breast cancer samples
- 46 healthy controls
- 92 total samples
- Blood plasma
- Whole-genome sequencing
- HiSeq 2000
- hg38 fragment coordinates

Raw fragment-level `.bgz` files are not included in the repository because of their large size.

## Workflow

1. Download fragment-level data from FinaleDB.
2. Process raw fragment files using `src/process_fragments.py`.
3. Apply MAPQ >= 30 filtering.
4. Restrict regional analysis to autosomes.
5. Generate global fragment-length distributions.
6. Calculate short and long fragment counts in 5-Mb genomic bins.
7. Construct global fragmentomic features.
8. Log-transform and center regional short/long ratios within each sample.
9. Compare global, regional, and combined classification models.
10. Simulate reduced sequencing depth using binomial thinning.
11. Measure regional-profile stability relative to full depth.
12. Evaluate breast cancer classification across simulated fragment counts.

## Main Results

Full-depth classification:

- Global features: ROC-AUC 0.646 ± 0.114
- Regional features: ROC-AUC 0.917 ± 0.056
- Combined features: ROC-AUC 0.921 ± 0.057

Replicate-wise simulated-depth classification:

- 20M fragments: ROC-AUC 0.871 ± 0.030
- 10M fragments: ROC-AUC 0.822 ± 0.037
- 5M fragments: ROC-AUC 0.803 ± 0.054
- 2M fragments: ROC-AUC 0.655 ± 0.079
- 1M fragments: ROC-AUC 0.619 ± 0.072

Regional-profile stability progressively decreased as sequencing depth was reduced.

## Repository Structure

- `data/processed/` — processed cohort-level data
- `notebooks/` — downstream analysis notebook
- `src/` — raw fragment preprocessing script
- `results/figures/` — final figures
- `results/tables/` — final result tables
- `report/` — final project report
- `environment.yml` — reproducible Conda environment

## Environment Setup

Create the environment with:

    conda env create -f environment.yml
    conda activate cfdna-fragmentomics

## Raw Fragment Processing

Place the downloaded hg38 FinaleDB fragment files into breast cancer and healthy input folders, then run:

    python src/process_fragments.py

The preprocessing script generates QC summaries, global fragment-length counts, and regional short/long fragment counts.

## Downstream Analysis

Open:

`notebooks/cfDNA_Fragmentomics_Final_Analysis.ipynb`

The notebook performs QC, feature construction, statistical testing, PCA, classification, simulated sequencing-depth reduction, feature-stability analysis, and final visualization.

## Reproducibility

Fixed random seeds were used for simulated downsampling and model evaluation. Raw sequencing-fragment files are not stored in GitHub because of their size. Processed cohort-level tables are provided where practical to support reproduction of the downstream analysis.

## Key Software

- Python 3.12
- NumPy
- pandas
- SciPy
- scikit-learn
- matplotlib
- statsmodels

## Author

Nelly Abualkheir

## Course

Week 8 Final Project
