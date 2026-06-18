# Figure 4 — Illumina EPIC Array Coverage and ART Data Analysis

This directory contains the code to reproduce Figure 4, which evaluates how well commercial Illumina EPIC arrays capture CoRSIVs and applies CoRSIVs to assisted reproductive technology (ART) data:

- Matching control genomic regions to CoRSIVs by CpG count and GC content
- Overlap between CoRSIV loci and Illumina EPIC v1/v2 probes
- Analysis of ART (IVF/ICSI) cohort methylation data at CoRSIV loci

---

## Subdirectories

### `01_control_region_matching/`

Generates a matched control set of genomic regions with similar CpG density and GC content to CoRSIVs for use as a null comparator throughout the paper.

**Inputs:**
- CoRSIV BED file (from Figure 1)
- hg38 CpG annotation

**Run:**
```bash
python match_control_probes.py
python plot_matching_qc.py        # QC plots verifying the match quality
python plot_nearest_2d_matching.py
```

---

### `02_EPIC_probe_overlap/`

Calculates what fraction of CoRSIVs contain at least one Illumina EPIC probe.

**Inputs:**
- CoRSIV BED file
- Illumina EPIC manifest (place in `data/raw/`)

**Run:**
```bash
bash make_EPIC_CoRSIV_table.sh
python make_EPIC_CoRSIV_table.py
```

---

### `03_ART_data_analysis/`

Applies CoRSIV BED files to filter EPIC array data from ART cohorts, then plots methylation differences.

**Inputs:**
- EPIC beta-value matrix from ART cohort (place in `data/raw/`)
- CoRSIV BED file

**Run:**
```bash
python filter_by_ID_REF.py
python remove_detection_pval_cols.py
```

---

## Dependencies

- Python: `pandas`, `numpy`, `matplotlib`, `scipy`
- Shell: `bedtools ≥ 2.30`
