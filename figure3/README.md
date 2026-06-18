# Figure 3 — Cross-Tissue Systemic Interindividual Variation (SIV)

This directory contains the code to reproduce Figure 3, which validates that CoRSIVs exhibit correlated methylation across tissues within individuals:

- Pairwise tissue–tissue methylation correlation matrices
- Hierarchical clustering heatmap of individuals based on CoRSIV methylation
- mQTL analysis linking CoRSIV methylation to nearby genetic variants

---

## Subdirectories

### `01_tissue_correlation/`

Computes Pearson/Spearman correlation of methylation beta values between every tissue pair at CoRSIV loci and plots the correlation matrix.

**Inputs:**
- Multi-tissue methylation matrix (beta values, samples × CpGs) — place in `data/raw/`
- CoRSIV BED file (from Figure 1)

**Steps (run in order):**
```bash
python correlation_matrix_calc_step1.py
python correlation_heatmap_step2.py
```

---

### `02_SIV_clustering/`

Identifies CoRSIVs that show systemic interindividual variation (correlated across tissues) using region-wise tissue-pair statistics, then clusters individuals.

**Steps (run in order):**
```bash
python correlation_matrix_calc_step1.py        # shared with 01
jupyter nbconvert --to notebook --execute SIV_Analysis_methy_clustering_step3.ipynb
python calc_regionwise_tissuepair_step4.py
```

---

### `03_mQTL/`

Identifies methylation QTLs (mQTLs) at CoRSIV loci and plots the results.

**Inputs:**
- Genotype data (VCF/PLINK format) — place in `data/raw/`
- CoRSIV methylation matrix

**Run:**
```r
Rscript genetic_columns.R
Rscript plot_mqtl_fig.R
```

---

## Dependencies

- Python: `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`, `jupyter`
- R: `ggplot2`, `dplyr`, `MatrixEQTL` (for mQTL)
