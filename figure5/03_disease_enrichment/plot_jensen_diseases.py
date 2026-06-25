#!/usr/bin/env python3
"""
Combined Jensen DISEASES Enrichment Heatmap
===========================================
Creates a single compact heatmap panel for Jensen DISEASES terms.
Compares black vs white TF lists across cell lines.
Designed as a compact panel for a multi-panel figure.

Author: auto-generated
Date: 2026-03-12
"""

import os
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ── Configuration ──────────────────────────────────────────────────────────────

BASE_DIR = os.path.dirname(__file__)
RESULTS_DIR = os.path.join(BASE_DIR, "enrichr_results", "Jensen_DISEASES")
OUTPUT_DIR = os.path.join(BASE_DIR, "plots")
os.makedirs(OUTPUT_DIR, exist_ok=True)

CELL_LINES = ["HepG2", "K562", "MCF7", "WA01", "WA09"]
CATEGORIES = ["black", "white"]

# How many top terms to show (union of top black + top white, capped)
TOP_N_TERMS = 10
PVAL_THRESHOLD = 0.05

# Typography tuned for a compact 5.5 x 3 inch panel
TITLE_FS = 8
XTICK_FS = 6
TERM_FS = 5
CBAR_TICK_FS = 5
CBAR_LABEL_FS = 6
SUPTITLE_FS = 8


def format_term_two_lines(term, width=22):
    """Wrap term labels to at most two lines for compact center text."""
    wrapped = textwrap.wrap(term, width=width, break_long_words=False)
    if len(wrapped) <= 1:
        return term
    if len(wrapped) == 2:
        return "\n".join(wrapped)
    return f"{wrapped[0]}\n{' '.join(wrapped[1:])}"

# ── Data Loading ───────────────────────────────────────────────────────────────

def load_enrichr_results():
    """Load black and white Jensen DISEASES results across all cell lines."""
    all_dfs = []

    for category in CATEGORIES:
        for cell_line in CELL_LINES:
            filename = f"top300_TFs_{category}_{cell_line}__Jensen_DISEASES.tsv"
            filepath = os.path.join(RESULTS_DIR, filename)

            if not os.path.exists(filepath):
                print(f"  WARNING: {filename} not found, skipping.")
                continue

            df = pd.read_csv(filepath, sep="\t")
            df["Category"] = category
            df["CellLine"] = cell_line
            df["neg_log10_adjp"] = -np.log10(df["Adjusted_P-value"].clip(lower=1e-300))
            df["Term_clean"] = df["Term"].str.replace(r"\s+", " ", regex=True).str.strip()
            all_dfs.append(df)

    if not all_dfs:
        raise FileNotFoundError(
            f"No Jensen DISEASES TSV files were found in {RESULTS_DIR}"
        )

    combined = pd.concat(all_dfs, ignore_index=True)
    print(f"Loaded {len(combined)} total Jensen DISEASES enrichment entries")
    print(f"  Categories: {combined['Category'].unique()}")
    print(f"  Cell lines: {combined['CellLine'].unique()}")
    return combined


# ── Compact Jensen DISEASES Heatmap ───────────────────────────────────────────

def plot_jensen_heatmap(df):
    """
    Create a compact heatmap:
    Black (left) | disease labels (center) | White (right).
    """
    # Select top terms from each category, then union while preserving order.
    black_df = df[df["Category"] == "black"]
    top_black = (
        black_df.groupby("Term_clean")["neg_log10_adjp"]
        .max()
        .nlargest(TOP_N_TERMS)
        .index.tolist()
    )

    white_df = df[df["Category"] == "white"]
    top_white = (
        white_df.groupby("Term_clean")["neg_log10_adjp"]
        .max()
        .nlargest(TOP_N_TERMS)
        .index.tolist()
    )

    top_terms = list(dict.fromkeys(top_black + top_white))[:TOP_N_TERMS]
    print(f"  Total terms in panel: {len(top_terms)}")

    filtered = df[df["Term_clean"].isin(top_terms)]

    black_pivot = (
        filtered[filtered["Category"] == "black"]
        .pivot_table(
            index="Term_clean", columns="CellLine", values="neg_log10_adjp", aggfunc="first"
        )
        .reindex(index=top_terms, columns=CELL_LINES)
        .fillna(0)
    )

    white_pivot = (
        filtered[filtered["Category"] == "white"]
        .pivot_table(
            index="Term_clean", columns="CellLine", values="neg_log10_adjp", aggfunc="first"
        )
        .reindex(index=top_terms, columns=CELL_LINES)
        .fillna(0)
    )

    vmax = max(black_pivot.max().max(), white_pivot.max().max())
    vmax = max(vmax, 2)

    fig = plt.figure(figsize=(5.5, 3))
    gs = fig.add_gridspec(1, 4, width_ratios=[3, 7, 3, 0.3], wspace=0.02)
    ax1 = fig.add_subplot(gs[0, 0])
    ax_labels = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2])
    cax = fig.add_subplot(gs[0, 3])

    cmap = sns.color_palette("YlOrRd", as_cmap=True)

    sns.heatmap(
        black_pivot, ax=ax1, cmap=cmap, vmin=0, vmax=vmax,
        cbar=False, linewidths=0.8, linecolor="white",
        annot=False,
        xticklabels=True, yticklabels=False
    )
    ax1.set_title("Black", fontsize=TITLE_FS, fontweight="bold", pad=4)
    ax1.set_ylabel("")
    ax1.set_xlabel("")
    ax1.tick_params(axis="x", labelsize=XTICK_FS, rotation=35)
    for label in ax1.get_xticklabels():
        label.set_horizontalalignment("right")

    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(0, len(top_terms))
    ax_labels.invert_yaxis()
    for i, term in enumerate(top_terms):
        display_term = format_term_two_lines(term)
        ax_labels.text(
            0.5,
            i + 0.5,
            display_term,
            ha="center",
            va="center",
            fontsize=TERM_FS,
            fontweight="medium",
            linespacing=1.0,
        )
    ax_labels.axis("off")

    sns.heatmap(
        white_pivot, ax=ax2, cmap=cmap, vmin=0, vmax=vmax,
        cbar=True, cbar_ax=cax,
        cbar_kws={"label": "-log₁₀(Adj. P)"},
        linewidths=0.8, linecolor="white",
        annot=False,
        xticklabels=True, yticklabels=False
    )
    ax2.set_title("White", fontsize=TITLE_FS, fontweight="bold", pad=4)
    ax2.set_ylabel("")
    ax2.set_xlabel("")
    ax2.tick_params(axis="x", labelsize=XTICK_FS, rotation=35)
    for label in ax2.get_xticklabels():
        label.set_horizontalalignment("right")

    cax.tick_params(labelsize=CBAR_TICK_FS)
    cax.yaxis.label.set_size(CBAR_LABEL_FS)

    sig_line = -np.log10(PVAL_THRESHOLD)
    fig.text(
        0.5,
        0.02,
        f"Significance guide: -log₁₀(Adj. P) > {sig_line:.2f} (Adj. P < {PVAL_THRESHOLD})",
        ha="center",
        fontsize=CBAR_TICK_FS,
        color="#666666",
    )

    fig.suptitle(
        "Jensen DISEASES Enrichment: Black vs White TF Lists",
        fontsize=SUPTITLE_FS,
        fontweight="bold",
        y=0.99,
    )

    plt.subplots_adjust(left=0.03, right=0.93, top=0.84, bottom=0.24)

    outpath_png = os.path.join(OUTPUT_DIR, "jensen_diseases_combined_heatmap_black_vs_white.png")
    outpath_pdf = os.path.join(OUTPUT_DIR, "jensen_diseases_combined_heatmap_black_vs_white.pdf")
    fig.savefig(outpath_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(outpath_pdf, bbox_inches="tight", facecolor="white")
    print(f"\nCombined Jensen DISEASES heatmap saved: {outpath_png}")
    print(f"Combined Jensen DISEASES heatmap saved: {outpath_pdf}")
    plt.close(fig)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("Jensen DISEASES Enrichment — Black vs White Comparison")
    print("=" * 60)

    df = load_enrichr_results()

    print(f"\n--- Generating Combined Jensen Heatmap (top {TOP_N_TERMS}) ---")
    plot_jensen_heatmap(df)

    print("\n✅ Done!")


if __name__ == "__main__":
    main()
