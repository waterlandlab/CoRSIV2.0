# Figure 2 — Genomic Distribution and Population-Genetic Properties

This directory contains the code to reproduce Figure 2, which characterises the genomic landscape of CoRSIVs:

- Circular genome-wide density plot (circos-style)
- SNP enrichment counts around CoRSIV loci vs. control regions
- Subtelomeric enrichment bar plot

---

## Subdirectories

### `01_genome_distribution/`

Circular chromosome density plot showing the distribution of CoRSIVs across the genome.

**Inputs:** CoRSIV BED file (from Figure 1 `01_CoRSIV_screening_and_bed_generation/`)

**Run:**
```r
Rscript circlize_density_plot.R
```

**Dependencies:** R packages `circlize`, `GenomicRanges`

---

### `02_SNP_enrichment/`

Counts SNPs from a VCF (e.g., 1000 Genomes) within windows around CoRSIVs and matched control regions, then plots the comparison.

**Inputs:**
- CoRSIV BED file
- Control region BED file
- 1000 Genomes VCF (place in `data/raw/`)

**Run:**
```r
Rscript make_comparison_plot.R
```

---

### `03_subtelomeric/`

Tests whether CoRSIVs are enriched in subtelomeric regions and generates the bar plot.

**Inputs:** CoRSIV BED file, hg38 chromosome arm coordinates

**Run:**
```r
Rscript sub_telomeric_bar_plot.R
```

---

## Dependencies

- R: `ggplot2`, `dplyr`, `circlize`, `GenomicRanges`
- BED files from Figure 1 must be generated before running these scripts
