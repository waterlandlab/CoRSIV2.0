# Figure 1 — CoRSIV Identification and Basic Characterisation

This directory contains the code to reproduce Figure 1, which covers:

- Filtering CoRSIV screening outputs by CpG count and interindividual range thresholds
- Generating BED files for Black, White, and multiethnic cohorts
- Venn diagrams showing overlap between population-specific CoRSIV sets
- Methylation heatmap and PCA across samples
- CpG density analysis comparing CoRSIVs to the rest of the genome

---

## Subdirectories

### `01_CoRSIV_screening_and_bed_generation/`

Filters raw screening CSVs and produces the BED files used by all downstream analyses.

**Inputs** (`data/raw/`):
- `black_CoRSIV_screen_output.csv`
- `white_CoRSIV_screen_output.csv`
- `multiracial_CoRSIV_screen_output.csv`

**Outputs** (`data/processed/beds_by_condition_*/`):
- `CpGge{3|4|5}_Rangegt{0|10|20}.bed` — one BED file per filtering condition (3 × 3 = 9 files per cohort)

**Run:**
```bash
jupyter nbconvert --to notebook --execute plots_and_generate_bed_files_2.ipynb --output plots_and_generate_bed_files_2_executed.ipynb
```

**Key filtering parameters** (edit at the top of the notebook):
| Parameter | Values tested | Default for figures |
|---|---|---|
| CpG Count Sum | ≥ 3, ≥ 4, ≥ 5 | ≥ 5 |
| Interindividual range | > 0%, > 10%, > 20% | > 20% |

---

### `02_venn_diagram/`

Compares CoRSIV sets across populations using Venn diagrams and generates a unified consensus track.

**Inputs:** BED files from `01_CoRSIV_screening_and_bed_generation/`

**Run:**
```bash
bash make_unified_track.sh
Rscript venn.R
```

---

### `03_heatmap_PCA/`

Methylation heatmap with sample annotations and PCA plot.

**Inputs:** Methylation beta-value matrix (place in `data/raw/`)

**Run:**
```bash
Rscript create_annotated_heatmap.R
python PCA_plot_CORSIV_samples.py
```

---

### `04_CpG_density/`

Calculates and plots CpG density within CoRSIVs vs. genome background.

**Requires:** `bedtools`, hg38 reference FASTA

**Run:**
```bash
bash download_hg38.sh          # downloads reference if not present
bash run_bedtools_nuc.sh       # computes nucleotide composition per region
python plot_cpg_density_with_islands.py
```

---

## Dependencies

- Python: `pandas`, `matplotlib`, `numpy`, `jupyter`
- R: `ggplot2`, `dplyr`, `VennDiagram`
- Shell: `bedtools ≥ 2.30`, `wget`
