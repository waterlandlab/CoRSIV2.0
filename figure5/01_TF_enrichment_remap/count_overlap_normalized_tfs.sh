#!/usr/bin/env bash
set -euo pipefail

outdir="overlap_results_norm"
mkdir -p "${outdir}"
export MPLCONFIGDIR="${outdir}/.mplconfig"
mkdir -p "${MPLCONFIGDIR}"

counts_tsv="${outdir}/tf_overlap_normalized_counts.tsv"
echo -e "query_file\tremap_file\tcell_line\ttf\tquery_count\tnormalized_count" > "${counts_tsv}"

shopt -s nullglob
query_files=(query/*.bed)
remap_files=(remap2022_all_macs2_hg38_v1_0_*.bed)

if [[ ${#query_files[@]} -eq 0 ]]; then
  echo "No query BED files found in query/" >&2
  exit 1
fi

if [[ ${#remap_files[@]} -eq 0 ]]; then
  echo "No remap BED files found in current directory" >&2
  exit 1
fi

for q in "${query_files[@]}"; do
  qbase=$(basename "${q}")
  if [[ "${qbase}" == "sorted.matched.controls.bed" || "${qbase}" == "hg38.chrom.sizes.txt" ]]; then
    continue
  fi
  
  # Calculate total regions in this query file (excluding empty lines if any)
  total_regions=$(awk 'NF' "${q}" | wc -l | tr -d ' ')
  
  if [[ ${total_regions} -eq 0 ]]; then
    echo "Warning: ${q} is empty, skipping." >&2
    continue
  fi

  for r in "${remap_files[@]}"; do
    rbase=$(basename "${r}")
    cell=$(echo "${rbase}" | sed -E 's/^remap2022_all_macs2_hg38_v1_0_//' | sed -E 's/\.bed$//')

    # Count unique query intervals overlapping each TF (not remap peak counts).
    bedtools intersect -a "${q}" -b "${r}" -wa -wb \
      | awk -v q="${qbase}" -v r="${rbase}" -v c="${cell}" -v N="${total_regions}" 'BEGIN{FS=OFS="\t"}{
          tf_field=$(NF-5);
          n=split(tf_field, a, ".");
          tf=(n>=2 ? a[2] : tf_field);
          qkey=$1 ":" $2 ":" $3;
          seen[tf, qkey]=1;
        } END{
          for (k in seen) {
            split(k, parts, SUBSEP);
            tf=parts[1];
            counts[tf]++;
          }
          for (tf in counts) {
            norm_count = (counts[tf] / N) * 1000;
            print q, r, c, tf, counts[tf], norm_count;
          }
        }' \
      >> "${counts_tsv}"
  done
done

python - <<'PY'
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

df = pd.read_csv("overlap_results_norm/tf_overlap_normalized_counts.tsv", sep="\t")

if df.empty:
    raise SystemExit("No overlaps detected; plot not created.")

# Keep top 10 TFs in each query/remap pair based on normalized_count.
top = (
    df.sort_values(["query_file", "remap_file", "normalized_count"], ascending=[True, True, False])
      .groupby(["query_file", "remap_file"], as_index=False)
      .head(10)
      .copy()
)

top["panel"] = top["query_file"] + " | " + top["cell_line"]

# Sort TF labels within each panel by count.
top = top.sort_values(["panel", "normalized_count"], ascending=[True, False])

sns.set_theme(style="whitegrid")
g = sns.catplot(
    data=top,
    x="normalized_count",
    y="tf",
    col="panel",
    col_wrap=5,
    kind="bar",
    sharex=False,
    sharey=False,
    height=3.6,
    aspect=1.15,
)
g.set_axis_labels("Normalized Overlaps (per 1000 regions)", "TF")
g.set_titles("{col_name}")
for ax in g.axes.flatten():
    ax.tick_params(axis="y", labelsize=8)

plt.tight_layout()
plt.savefig("overlap_results_norm/top10_tfs_overlap_norm_counts.png", dpi=300)

# Build per-cell-line TF x query-file heatmaps with row-wise min-max normalization.
heatmap_dir = Path("overlap_results_norm/heatmaps")
heatmap_dir.mkdir(parents=True, exist_ok=True)

for cell_line, sub in df.groupby("cell_line"):
    mat = (
        sub.groupby(["tf", "query_file"], as_index=False)["normalized_count"]
           .sum()
           .pivot(index="tf", columns="query_file", values="normalized_count")
           .fillna(0)
    )
    if mat.empty:
        continue

    # Keep top 50 TFs by total normalized counts in this cell line.
    top_tfs = mat.sum(axis=1).sort_values(ascending=False).head(50).index
    mat = mat.loc[top_tfs]

    # Normalize each TF row to [0, 1] across query files.
    row_min = mat.min(axis=1)
    row_max = mat.max(axis=1)
    denom = (row_max - row_min).replace(0, 1)
    norm = mat.sub(row_min, axis=0).div(denom, axis=0)

    # Sort rows by normalized mean and then max for a stable, interpretable layout.
    order = norm.mean(axis=1).sort_values(ascending=False).index
    norm = norm.loc[order]
    mat = mat.loc[order]

    # Save both raw and normalized matrices.
    cell_slug = cell_line.replace("/", "_")
    mat.to_csv(heatmap_dir / f"tf_query_matrix_raw_norm_{cell_slug}.tsv", sep="\t")
    norm.to_csv(heatmap_dir / f"tf_query_matrix_norm01_{cell_slug}.tsv", sep="\t")

    # Dynamic height: readable without creating an enormous image.
    h = max(4, min(40, 0.16 * len(norm) + 1.5))
    w = max(6, 1.2 * len(norm.columns) + 2.5)
    plt.figure(figsize=(w, h))
    sns.heatmap(
        norm,
        cmap="mako",
        vmin=0,
        vmax=1,
        cbar_kws={"label": "Row-normalized query count (0-1)"},
        linewidths=0.0,
    )
    plt.title(f"TF vs Query File Heatmap (Normalized 0-1): {cell_line}")
    plt.xlabel("Query file")
    plt.ylabel("TF")
    plt.tight_layout()
    plt.savefig(heatmap_dir / f"tf_query_heatmap_norm01_{cell_slug}.png", dpi=300)
    plt.close()
PY

echo "Wrote:"
echo "  ${counts_tsv}"
echo "  ${outdir}/top10_tfs_overlap_norm_counts.png"
echo "  ${outdir}/heatmaps/tf_query_heatmap_norm01_<cell_line>.png"
