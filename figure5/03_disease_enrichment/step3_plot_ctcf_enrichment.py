#!/usr/bin/env python3
"""
Plot CTCF dominance among highly significant overlap enrichments.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


INPUT_CSV_NAME = "all_results_overlap_heatmap_data_wide.csv"
OUTPUT_THRESHOLD_FIG = Path("ctcf_enrichment_threshold_share.png")
OUTPUT_THRESHOLD_PDF = Path("ctcf_enrichment_threshold_share.pdf")
OUTPUT_SUMMARY = Path("ctcf_enrichment_threshold_summary.csv")

PVALUE_SCORE_COL = "max_neglog10_q"  # larger means more significant
TF_COL = "TF"
EXPERIMENT_COL = "experiment_id"
SIGNIFICANT_THRESHOLD = 20
TOP_TF_THRESHOLDS = [10, 20, 30, 40, 50]
OTHER_TFS_TO_PLOT = ["RAD21", "STAG1", "SMC3", "STAG2", "SMC1A"]

TITLE_FONTSIZE = 20
LABEL_FONTSIZE = 18
TICK_FONTSIZE = 15
LEGEND_FONTSIZE = 14
ANNOTATION_FONTSIZE = 13
SUPTITLE_FONTSIZE = 19


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    input_candidates = [
        Path(INPUT_CSV_NAME),
        script_dir / INPUT_CSV_NAME,
        script_dir / "remapenrich_results_2022_with_universe" / INPUT_CSV_NAME,
    ]
    input_csv = next((p for p in input_candidates if p.exists()), None)
    if input_csv is None:
        raise FileNotFoundError(f"Could not find {INPUT_CSV_NAME} in expected locations.")

    df = pd.read_csv(input_csv)
    if TF_COL not in df.columns:
        if "category" not in df.columns:
            raise KeyError(f"Missing both '{TF_COL}' and 'category' columns in input CSV.")
        parts = df["category"].astype(str).str.split(".")
        df[TF_COL] = parts.str[1].where(parts.str.len() > 1, parts.str[0])

    if EXPERIMENT_COL not in df.columns:
        if "category" in df.columns:
            df[EXPERIMENT_COL] = df["category"].astype(str).str.split(".").str[0].str.strip()
        else:
            raise KeyError(f"Missing both '{EXPERIMENT_COL}' and 'category' columns in input CSV.")

    df = df[[EXPERIMENT_COL, TF_COL, PVALUE_SCORE_COL]].dropna().copy()
    df[TF_COL] = df[TF_COL].astype(str).str.strip()
    df[EXPERIMENT_COL] = df[EXPERIMENT_COL].astype(str).str.strip()
    df["is_ctcf"] = df[TF_COL].str.upper().eq("CTCF")

    # Panel A: top TFs within highly significant overlaps, counted as unique experiments.
    sig = df[df[PVALUE_SCORE_COL] >= SIGNIFICANT_THRESHOLD].copy()
    top_counts = (
        sig.groupby(TF_COL)[EXPERIMENT_COL]
        .nunique()
        .sort_values(ascending=False)
        .head(10)
    )
    colors = ["#d62728" if tf.upper() == "CTCF" else "#7f8c8d" for tf in top_counts.index]

    # Panel B: CTCF and selected other TF fractions as threshold tightens, using unique experiments.
    thresholds = [10, 20, 30, 40, 50, 60, 70, 80]
    rows = []
    for thr in thresholds:
        sub = df[df[PVALUE_SCORE_COL] >= thr]
        n = int(sub[EXPERIMENT_COL].nunique())
        ctcf_n = int(sub.loc[sub["is_ctcf"], EXPERIMENT_COL].nunique())
        frac = (ctcf_n / n) if n else 0.0
        row = {
            "threshold": thr,
            "n_unique_experiments": n,
            "ctcf_unique_experiments": ctcf_n,
            "ctcf_fraction": frac,
        }
        for tf in OTHER_TFS_TO_PLOT:
            tf_count = int(
                sub.loc[sub[TF_COL].str.upper().eq(tf.upper()), EXPERIMENT_COL].nunique()
            )
            row[f"{tf}_unique_experiments"] = tf_count
            row[f"{tf}_fraction"] = (tf_count / n) if n else 0.0
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(OUTPUT_SUMMARY, index=False)

    # Figure 1: Top TFs bar plots for selected significance cutoffs.
    top_tf_outputs = []
    for top_tf_thr in TOP_TF_THRESHOLDS:
        top_sig = df[df[PVALUE_SCORE_COL] >= top_tf_thr].copy()
        top_counts = (
            top_sig.groupby(TF_COL)[EXPERIMENT_COL]
            .nunique()
            .sort_values(ascending=False)
            .head(10)
        )
        fig1, ax1 = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
        ax1.set_title(f"Top TFs at max_neglog10_q >= {top_tf_thr}", fontsize=TITLE_FONTSIZE)
        ax1.set_ylabel("Count of unique experiments", fontsize=LABEL_FONTSIZE)
        if top_counts.empty:
            ax1.text(
                0.5,
                0.5,
                "No experiments at this threshold",
                ha="center",
                va="center",
                transform=ax1.transAxes,
                fontsize=LABEL_FONTSIZE,
                color="#2f3640",
            )
            ax1.set_xticks([])
            ax1.set_yticks([])
        else:
            colors = ["#d62728" if tf.upper() == "CTCF" else "#7f8c8d" for tf in top_counts.index]
            ax1.bar(top_counts.index, top_counts.values, color=colors, edgecolor="black", linewidth=0.4)
            ax1.tick_params(axis="x", rotation=55, labelsize=TICK_FONTSIZE)
            ax1.tick_params(axis="y", labelsize=TICK_FONTSIZE)
            max_count = int(top_counts.max())
            for i, val in enumerate(top_counts.values):
                ax1.text(
                    i,
                    val + max(1, 0.015 * max_count),
                    str(int(val)),
                    ha="center",
                    va="bottom",
                    fontsize=ANNOTATION_FONTSIZE,
                )

        out_png = Path(f"ctcf_enrichment_top_tfs_thr{top_tf_thr}.png")
        out_pdf = Path(f"ctcf_enrichment_top_tfs_thr{top_tf_thr}.pdf")
        fig1.savefig(out_png, dpi=300, bbox_inches="tight")
        fig1.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig1)
        top_tf_outputs.extend([out_png, out_pdf])

    # Figure 2: Threshold share line plot.
    fig2, ax2 = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    ax2.plot(
        summary["threshold"],
        summary["ctcf_fraction"] * 100,
        marker="o",
        color="#d62728",
        linewidth=2.4,
        label="CTCF",
        zorder=3,
    )
    for tf in OTHER_TFS_TO_PLOT:
        frac_col = f"{tf}_fraction"
        ax2.plot(
            summary["threshold"],
            summary[frac_col] * 100,
            marker="o",
            color="#9aa0a6",
            linewidth=1.2,
            alpha=0.9,
            label=tf,
            zorder=2,
        )

    ax2.set_title("TF share among unique experiments above threshold", fontsize=TITLE_FONTSIZE)
    ax2.set_xlabel("max_neglog10_q threshold", fontsize=LABEL_FONTSIZE)
    ax2.set_ylabel("Fraction of unique experiments (%)", fontsize=LABEL_FONTSIZE)
    ax2.set_xlim(0, 30)
    ax2.set_xticks([10, 20, 30])
    ax2.tick_params(axis="both", labelsize=TICK_FONTSIZE)
    ax2.set_ylim(0, 100)
    ax2.grid(axis="y", linestyle="--", alpha=0.4)
    ax2.legend(frameon=False, fontsize=LEGEND_FONTSIZE, ncols=2, loc="upper left")

    # Annotate sample sizes for each threshold.
    for _, r in summary.iterrows():
        ax2.annotate(
            f"n={int(r['n_unique_experiments'])}",
            (r["threshold"], r["ctcf_fraction"] * 100),
            textcoords="offset points",
            xytext=(0, 8),
            ha="center",
            fontsize=ANNOTATION_FONTSIZE,
            color="#2f3640",
        )

    overall_ctcf = int(sig.loc[sig["is_ctcf"], EXPERIMENT_COL].nunique())
    overall_n = int(sig[EXPERIMENT_COL].nunique())
    overall_pct = (100 * overall_ctcf / overall_n) if overall_n else 0.0
    fig2.suptitle(
        f"CTCF is the dominant TF in highly significant overlaps "
        f"({overall_ctcf}/{overall_n}, {overall_pct:.1f}% at threshold >= {SIGNIFICANT_THRESHOLD})",
        fontsize=SUPTITLE_FONTSIZE,
        y=1.01,
    )

    fig2.savefig(OUTPUT_THRESHOLD_FIG, dpi=300, bbox_inches="tight")
    fig2.savefig(OUTPUT_THRESHOLD_PDF, bbox_inches="tight")
    plt.close(fig2)
    for output_path in top_tf_outputs:
        print(f"Wrote figure: {output_path}")
    print(f"Wrote figure: {OUTPUT_THRESHOLD_FIG}")
    print(f"Wrote figure: {OUTPUT_THRESHOLD_PDF}")
    print(f"Wrote summary: {OUTPUT_SUMMARY}")


if __name__ == "__main__":
    main()
