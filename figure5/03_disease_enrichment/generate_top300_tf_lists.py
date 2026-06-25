#!/usr/bin/env python3
"""
Generate Top-300 TF Gene Lists
================================
Reads tf_overlap_normalized_counts.tsv and writes one .txt file per
(query_type × cell_line) combination into the top300_tf_lists/ folder.

Each file contains up to 300 TF gene symbols ranked by normalized_count
(descending), one symbol per line.

Output naming convention:
    top300_TFs_<query>_<cell_line>.txt
    e.g. top300_TFs_black_HepG2.txt
         top300_TFs_white_K562.txt
         top300_TFs_random_MCF7.txt
"""

import os
import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_TSV  = os.path.join(SCRIPT_DIR, "tf_overlap_normalized_counts.tsv")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "top300_tf_lists")
TOP_N      = 300

# ── Query-file → short label mapping ──────────────────────────────────────────
QUERY_LABELS = {
    "Black_CpGge5_Rangegt20.bed": "black",
    "White_CpGge5_Rangegt20.bed": "white",
    "random_background.bed":      "random",
}

# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    print(f"Reading {INPUT_TSV} ...")
    df = pd.read_csv(INPUT_TSV, sep="\t")

    # Sanity-check expected columns
    required = {"query_file", "cell_line", "tf", "normalized_count"}
    missing  = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in TSV: {missing}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    cell_lines   = sorted(df["cell_line"].unique())
    query_files  = sorted(df["query_file"].unique())

    total_files = 0

    for query_file in query_files:
        # Map the raw filename to a clean label; skip unknowns
        query_label = QUERY_LABELS.get(query_file)
        if query_label is None:
            print(f"  WARNING: unrecognised query_file '{query_file}' — skipping.")
            continue

        for cell_line in cell_lines:
            subset = df[
                (df["query_file"] == query_file) &
                (df["cell_line"]  == cell_line)
            ].copy()

            if subset.empty:
                print(f"  No data for {query_label} / {cell_line} — skipping.")
                continue

            # Rank by normalized_count descending, take top N
            top = (
                subset
                .sort_values("normalized_count", ascending=False)
                .drop_duplicates(subset="tf")
                .head(TOP_N)
            )

            out_name = f"top{TOP_N}_TFs_{query_label}_{cell_line}.txt"
            out_path = os.path.join(OUTPUT_DIR, out_name)

            with open(out_path, "w") as fh:
                fh.write("\n".join(top["tf"].tolist()) + "\n")

            print(f"  Wrote {len(top):3d} TFs → {out_name}")
            total_files += 1

    print(f"\nDone! {total_files} file(s) saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
