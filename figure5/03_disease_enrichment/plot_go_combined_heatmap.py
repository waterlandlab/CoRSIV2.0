#!/usr/bin/env python3
"""
Combined GO Enrichment Heatmap
===============================
Creates a single heatmap combining top terms from:
  - GO Biological Process 2023
  - GO Molecular Function 2023
  - GO Cellular Component 2023

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
import matplotlib.patches as mpatches
import seaborn as sns

# ── Configuration ──────────────────────────────────────────────────────────────

BASE_DIR = os.path.dirname(__file__)
RESULTS_DIR = os.path.join(BASE_DIR, "enrichr_results")
OUTPUT_DIR = os.path.join(BASE_DIR, "plots")
os.makedirs(OUTPUT_DIR, exist_ok=True)

CELL_LINES = ["HepG2", "K562", "MCF7", "WA01", "WA09"]
CATEGORIES = ["black", "white"]

GO_LIBRARIES = {
    "GO_Biological_Process_2023": "Biological Process",
    "GO_Molecular_Function_2023": "Molecular Function",
    "GO_Cellular_Component_2023": "Cellular Component",
}

# How many top terms per GO library
TOP_N_PER_LIB = 3
PVAL_THRESHOLD = 0.05

# Typography tuned for a compact 5.5 x 3 inch panel
TITLE_FS = 8
XTICK_FS = 6
TERM_FS = 6.5
CBAR_TICK_FS = 5
CBAR_LABEL_FS = 6
LEGEND_FS = 5
SUPTITLE_FS = 8


def format_term_two_lines(term, width=22):
    """Wrap term labels to at most two lines for compact center text."""
    wrapped = textwrap.wrap(term, width=width, break_long_words=False)
    if len(wrapped) <= 1:
        return term
    if len(wrapped) == 2:
        return "\n".join(wrapped)
    # Collapse extra wrapped lines into two lines.
    return f"{wrapped[0]}\n{' '.join(wrapped[1:])}"


# ── Data Loading ───────────────────────────────────────────────────────────────

def load_go_results():
    """Load black and white results for all three GO libraries."""
    all_dfs = []

    for lib_key, lib_label in GO_LIBRARIES.items():
        lib_dir = os.path.join(RESULTS_DIR, lib_key)
        if not os.path.isdir(lib_dir):
            print(f"  WARNING: directory {lib_key} not found, skipping.")
            continue

        for category in CATEGORIES:
            for cell_line in CELL_LINES:
                filename = f"top300_TFs_{category}_{cell_line}__{lib_key}.tsv"
                filepath = os.path.join(lib_dir, filename)

                if not os.path.exists(filepath):
                    print(f"  WARNING: {filename} not found, skipping.")
                    continue

                df = pd.read_csv(filepath, sep="\t")
                df["Category"] = category
                df["CellLine"] = cell_line
                df["Library"] = lib_label
                df["neg_log10_adjp"] = -np.log10(
                    df["Adjusted_P-value"].clip(lower=1e-300)
                )
                # Clean up term names: remove GO ID in parentheses
                df["Term_clean"] = df["Term"].str.replace(
                    r"\s*\(GO:\d+\)\s*$", "", regex=True
                )
                all_dfs.append(df)

    combined = pd.concat(all_dfs, ignore_index=True)
    print(f"Loaded {len(combined)} total GO enrichment entries")
    print(f"  Libraries: {combined['Library'].unique()}")
    print(f"  Categories: {combined['Category'].unique()}")
    print(f"  Cell lines: {combined['CellLine'].unique()}")
    return combined


# ── Compact Combined GO Heatmap ───────────────────────────────────────────────

def plot_combined_go_heatmap(df):
    """
    Create a compact combined heatmap with top N terms from each GO library.
    Layout: Black (left) | term labels (center) | White (right)
    Rows are grouped by GO library with visual separators.
    """
    # Get top terms per library (union of black + white top terms)
    all_top_terms = []
    library_groups = []  # (library_label, [terms])

    for lib_label in GO_LIBRARIES.values():
        lib_df = df[df["Library"] == lib_label]
        if lib_df.empty:
            continue

        black_top = (
            lib_df[lib_df["Category"] == "black"]
            .groupby("Term_clean")["neg_log10_adjp"]
            .max()
            .nlargest(TOP_N_PER_LIB)
            .index.tolist()
        )
        white_top = (
            lib_df[lib_df["Category"] == "white"]
            .groupby("Term_clean")["neg_log10_adjp"]
            .max()
            .nlargest(TOP_N_PER_LIB)
            .index.tolist()
        )
        # Union, keep order, cap at TOP_N_PER_LIB
        top_terms = list(dict.fromkeys(black_top + white_top))[:TOP_N_PER_LIB]
        all_top_terms.extend(top_terms)
        library_groups.append((lib_label, top_terms))

    n_terms = len(all_top_terms)
    print(f"  Total terms across libraries: {n_terms}")

    # Build combined pivot tables
    black_rows = []
    white_rows = []

    for lib_label, terms in library_groups:
        lib_df = df[df["Library"] == lib_label]
        for term in terms:
            for cat, rows_list in [("black", black_rows), ("white", white_rows)]:
                cat_df = lib_df[
                    (lib_df["Category"] == cat) & (lib_df["Term_clean"] == term)
                ]
                row = {}
                for cl in CELL_LINES:
                    val = cat_df.loc[cat_df["CellLine"] == cl, "neg_log10_adjp"]
                    row[cl] = val.values[0] if len(val) > 0 else 0.0
                rows_list.append(row)

    black_pivot = pd.DataFrame(black_rows, index=all_top_terms, columns=CELL_LINES)
    white_pivot = pd.DataFrame(white_rows, index=all_top_terms, columns=CELL_LINES)

    # Shared color scale
    vmax = max(black_pivot.max().max(), white_pivot.max().max())
    vmax = max(vmax, 2)

    # ── Figure layout ──
    fig = plt.figure(figsize=(5.5, 3))
    gs = fig.add_gridspec(1, 4, width_ratios=[3.5, 5.0, 3.5, 0.3], wspace=0.005)
    ax1 = fig.add_subplot(gs[0, 0])
    ax_labels = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2])
    cax = fig.add_subplot(gs[0, 3])

    cmap = sns.color_palette("YlOrRd", as_cmap=True)

    # ── Black heatmap ──
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

    # ── White heatmap ──
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

    # ── Center disease labels with library group separators ──
    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(0, n_terms)
    ax_labels.invert_yaxis()

    # Library group colors for subtle background shading
    lib_colors = ["#f0f4ff", "#fff4f0", "#f0fff4"]

    row_idx = 0
    for grp_i, (lib_label, terms) in enumerate(library_groups):
        bg_color = lib_colors[grp_i % len(lib_colors)]
        # Background shading on label column
        ax_labels.axhspan(row_idx, row_idx + len(terms), color=bg_color, alpha=0.5)
        # Also shade the heatmap axes for visual grouping
        ax1.axhspan(row_idx, row_idx + len(terms), color=bg_color, alpha=0.15)
        ax2.axhspan(row_idx, row_idx + len(terms), color=bg_color, alpha=0.15)

        # Library name as a bracket/label on the side 
        mid_y = row_idx + len(terms) / 2
        
        for i, term in enumerate(terms):
            display_term = format_term_two_lines(term)
            ax_labels.text(0.5, row_idx + i + 0.5, display_term, ha="center", va="center",
                           fontsize=TERM_FS, fontweight="medium", linespacing=1.0)

        # Draw separator line between groups
        if grp_i < len(library_groups) - 1:
            sep_y = row_idx + len(terms)
            ax1.axhline(y=sep_y, color="#333333", linewidth=2)
            ax2.axhline(y=sep_y, color="#333333", linewidth=2)
            ax_labels.axhline(y=sep_y, color="#333333", linewidth=2)

        row_idx += len(terms)

    ax_labels.axis("off")

    # Colorbar font
    cax.tick_params(labelsize=CBAR_TICK_FS)
    cax.yaxis.label.set_size(CBAR_LABEL_FS)

    # ── Legend for GO library groups ──
    legend_patches = []
    for grp_i, (lib_label, _) in enumerate(library_groups):
        bg_color = lib_colors[grp_i % len(lib_colors)]
        legend_patches.append(
            mpatches.Patch(facecolor=bg_color, edgecolor="#999999",
                           label=lib_label, alpha=0.7)
        )
    fig.legend(handles=legend_patches, loc="lower center", ncol=len(library_groups),
               fontsize=LEGEND_FS, framealpha=0.9, edgecolor="#cccccc",
               bbox_to_anchor=(0.50, -0.08))

    fig.suptitle("GO Enrichment: Black vs White TF Lists",
                 fontsize=SUPTITLE_FS, fontweight="bold", y=0.99)

    plt.subplots_adjust(left=0.03, right=0.93, top=0.84, bottom=0.28)

    outpath_png = os.path.join(OUTPUT_DIR, "go_combined_heatmap_black_vs_white.png")
    outpath_pdf = os.path.join(OUTPUT_DIR, "go_combined_heatmap_black_vs_white.pdf")
    fig.savefig(outpath_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(outpath_pdf, bbox_inches="tight", facecolor="white")
    print(f"\nCombined GO heatmap saved: {outpath_png}")
    print(f"Combined GO heatmap saved: {outpath_pdf}")
    plt.close(fig)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("Combined GO Enrichment — Black vs White Comparison")
    print("=" * 60)

    df = load_go_results()

    print(f"\n--- Generating Combined GO Heatmap (top {TOP_N_PER_LIB} per library) ---")
    plot_combined_go_heatmap(df)

    print("\n✅ Done!")


if __name__ == "__main__":
    main()
