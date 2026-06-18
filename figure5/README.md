# Figure 5 — Functional Enrichment, Evolutionary Constraint, and Disease Associations

This directory contains the code to reproduce Figure 5, which interrogates the functional and evolutionary significance of CoRSIVs:

- Transcription factor (TF) binding enrichment using ReMap2022
- Tajima's D analysis from 1000 Genomes Project to assess evolutionary constraint
- Disease/phenotype enrichment via Enrichr
- TF–CoRSIV–gene network analysis

---

## Subdirectories

### `01_TF_enrichment_remap/`

Overlaps CoRSIV BED files with ReMap2022 ChIP-seq peak catalogue to identify enriched TFs, using matched control regions as the background.

**Inputs:**
- CoRSIV BED file (from Figure 1)
- Matched control BED file (from Figure 4 `01_control_region_matching/`)
- ReMap2022 merged peak BED (place in `data/raw/`)

**Steps (run in order):**
```bash
bash extract_corsiv_bed_regions.sh
bash generate_random_background.sh
bash count_overlap_normalized_tfs.sh
python extract_top300_tfs.py
python plot_scatter_normalized_matched_controls.py
Rscript generate_venn_diagrams.R
```

---

### `02_TajimaD_analysis/`

Calculates Tajima's D in windows overlapping CoRSIVs using 1000 Genomes phased VCFs and compares to genome-wide background.

**Inputs:**
- 1000 Genomes Project VCFs (CEU / YRI populations) — place in `data/raw/`
- CoRSIV BED file

**Run:**
```bash
bash run_1000g_tajd_pipeline.sh   # calls vcftools tajima_d per window
python plot_tajd_fixed.py
python plot_tajd_boxplot.py
```

**Requires:** `vcftools ≥ 0.1.16`, `tabix`

---

### `03_disease_enrichment/`

Tests whether genes near CoRSIVs are enriched in disease gene lists using Enrichr.

**Run:**
```python
python run_enrichr.py            # submits gene lists to Enrichr API
python plot_enrichr_scatter.py   # scatter plot of enrichment results
```

---

### `04_gene_network/`

Builds a TF–CoRSIV–gene regulatory network and runs cell-type–specific enrichment.

**Steps (run in order):**
```bash
python run_enrichr_combined.py
python run_enrichr_pooled.py
python run_enrichr_intersection.py
python plot_enrichr_pooled.py
python plot_enrichr_celltype.py
```

Jupyter notebooks for interactive exploration:
- `TF_CoRSIV_Analysis.ipynb`
- `TF_Gene_Association.ipynb`
- `Gene_Venn_Diagrams.ipynb`

---

## Dependencies

- Python: `pandas`, `numpy`, `matplotlib`, `scipy`, `requests` (Enrichr API)
- R: `ggplot2`, `dplyr`, `VennDiagram`, `remapenrich`
- Shell: `bedtools ≥ 2.30`, `vcftools ≥ 0.1.16`, `tabix`
