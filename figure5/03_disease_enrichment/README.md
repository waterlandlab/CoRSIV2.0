# Step 1: ReMapEnrich run with query + universe backgrounds

This step runs TF enrichment with `ReMapEnrich::enrichment()` using:
- Query BEDs: `query/Black_CpGge5_Rangegt20.bed`, `query/White_CpGge5_Rangegt20.bed`
- Universe BEDs: `background/sorted.matched.controls.bed`, `background/random_background.bed`
- Catalog BED: `remap2022_all_macs2_hg38_v1_0.bed`

Reference docs: [ReMapEnrich basic use](https://remap-cisreg.github.io/ReMapEnrich/vignettes/basic_use.html)

## What happened
The original script crashed on the `random_background.bed` runs with:

`NAs are not allowed in subscripted assignments`

This occurred inside `shuffleUniverse` during randomization.

## What we changed in Step 1
We updated `step1_run_remapenrich_with_universe_plots.R` to make preprocessing robust before calling `enrichment()`:

1. Remove invalid universe regions (NA/zero width).
2. Keep only query regions that can be placed in the universe (`query width <= max universe width`).
3. Restrict the universe to chromosomes shared with the query.
4. Record filtered counts in `run_summary.csv`.
5. Skip a run only if filtering leaves zero query or zero universe regions.

## Run command
```bash
Rscript step1_run_remapenrich_with_universe_plots.R
```

## Outputs
Output folder: `remapenrich_results_2022_with_universe/`

Per query-vs-universe pair:
- `*_enrichment.csv`
- `*_top20_by_qvalue.csv`
- `*_dotplot.pdf`

Combined summary:
- `run_summary.csv`

## Current status
Step 1 now completes successfully for all 4 query/universe combinations.

# Step 2: Aggregate Enrichment Results For Heatmap Input

## Manuscript-ready methods text (journal style)
We aggregated all four ReMapEnrich result tables (`*_enrichment.csv`) generated in Step 1 using [`step2_export_all_results_heatmap_data.R`](/Users/gunaseka/Documents/figure5_reanalysis/D_E/step2_export_all_results_heatmap_data.R). From each table, only `category` and multiple-testing-adjusted `q.value` were retained; entries with non-finite values or `q.value <= 0` were excluded. Enrichment significance was represented as `-log10(q.value)`. For each category within each run, duplicate entries were collapsed by taking the maximum `-log10(q.value)` (equivalently the minimum `q.value`) to obtain a single run-level score per category. The aggregated data were exported in long format (`all_results_overlap_heatmap_data_long.csv`) and then reshaped to wide format (`all_results_overlap_heatmap_data_wide.csv`), with one column per run and one row per category. A cross-run summary metric, `max_neglog10_q`, was computed as the maximum run-level score per category, and categories were ranked in descending order of this value for downstream prioritization.

# Step 3: CTCF-Centric Ranking And Threshold Analysis

## Manuscript-ready methods text (journal style)
CTCF-focused downstream analysis was performed using [`step3_plot_ctcf_enrichment.py`](/Users/gunaseka/Documents/figure5_reanalysis/D_E/step3_plot_ctcf_enrichment.py) on `all_results_overlap_heatmap_data_wide.csv`. TF labels were parsed from the second token of the `category` field and experiment identifiers from the first token (tokenization by `.`). To avoid double-counting repeated entries, all summary statistics were computed on unique experiment identifiers rather than raw overlap rows. Top-factor bar plots report the top 10 TFs by number of unique experiments above each significance cutoff (`max_neglog10_q >= 10, 20, 30, 40, 50`) and were exported as `ctcf_enrichment_top_tfs_thr*.{png,pdf}`. Threshold-dependent TF composition was then evaluated at `max_neglog10_q` cutoffs of 10, 20, 30, 40, 50, 60, 70, and 80 by calculating, at each cutoff, the fraction of unique experiments assigned to CTCF and to comparator factors (RAD21, STAG1, SMC3, STAG2, SMC1A), with denominator equal to all unique experiments passing that cutoff. These values were written to `ctcf_enrichment_threshold_summary.csv` and visualized in `ctcf_enrichment_threshold_share.{png,pdf}` (displayed with x-axis range 0-30 and ticks at 10, 20, and 30 for panel legibility). This step is descriptive and does not introduce additional inferential testing beyond the adjusted `q.value` framework carried forward from Step 1.
