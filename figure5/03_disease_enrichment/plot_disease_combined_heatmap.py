#!/usr/bin/env python3
"""
Combined disease enrichment heatmaps
===================================
Creates grouped heatmaps combining top disease terms from:
  - Jensen DISEASES
  - DisGeNET

Two matched panels are generated:
  1. HepG2 / K562 / MCF7 only
  2. WA01 / WA09 only

Each panel compares black vs white TF lists and groups rows by library.
"""

import os
import textwrap

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    import seaborn as sns
except ModuleNotFoundError:
    sns = None

# ── Configuration ──────────────────────────────────────────────────────────────

BASE_DIR = os.path.dirname(__file__)
RESULTS_DIR = os.path.join(BASE_DIR, "enrichr_results")
OUTPUT_DIR = os.path.join(BASE_DIR, "plots")
os.makedirs(OUTPUT_DIR, exist_ok=True)

CATEGORIES = ["black", "white"]

DISEASE_LIBRARIES = {
    "Jensen_DISEASES": "Jensen DISEASES",
    "DisGeNET": "DisGeNET",
}

PANELS = {
    "somatic": {
        "cell_lines": ["HepG2", "K562", "MCF7"],
        "display_labels": {"HepG2": "HepG2", "K562": "K562", "MCF7": "MCF7"},
        "title": "Disease Enrichment: Black vs White TF Lists",
        "subtitle": "HepG2, K562, and MCF7",
        "output_stub": "disease_combined_heatmap_black_vs_white_no_WA",
    },
    "stem": {
        "cell_lines": ["WA01", "WA09"],
        "display_labels": {"WA01": "WA1", "WA09": "WA9"},
        "title": "Disease Enrichment: Black vs White TF Lists",
        "subtitle": "WA1 and WA9",
        "output_stub": "disease_combined_heatmap_black_vs_white_WA_only",
    },
}

TOP_N_PER_LIB = 5
PVAL_THRESHOLD = 0.05

TITLE_FS = 8
XTICK_FS = 6
TERM_FS = 6
CBAR_TICK_FS = 5
CBAR_LABEL_FS = 6
LEGEND_FS = 5
SUPTITLE_FS = 8


def format_term_two_lines(term, width=24):
    """Keep disease term labels on a single line."""
    return term


def draw_heatmap(ax, data, cmap, vmin, vmax, cax=None, cbar_label=None):
    """Draw a heatmap with seaborn when available, else use matplotlib."""
    if sns is not None:
        return sns.heatmap(
            data,
            ax=ax,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            cbar=cax is not None,
            cbar_ax=cax,
            cbar_kws={"label": cbar_label} if cbar_label else None,
            linewidths=0.8,
            linecolor="white",
            annot=False,
            xticklabels=True,
            yticklabels=False,
        )

    image = ax.imshow(
        data.to_numpy(),
        aspect="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )
    ax.set_xticks(np.arange(len(data.columns)))
    ax.set_xticklabels(data.columns)
    ax.set_yticks([])
    ax.set_xticks(np.arange(-0.5, len(data.columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(data.index), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)

    if cax is not None:
        cbar = plt.colorbar(image, cax=cax)
        if cbar_label:
            cbar.set_label(cbar_label)

    return image


def load_disease_results():
    """Load black and white disease results for all configured libraries."""
    all_dfs = []

    for lib_key, lib_label in DISEASE_LIBRARIES.items():
        lib_dir = os.path.join(RESULTS_DIR, lib_key)
        if not os.path.isdir(lib_dir):
            print(f"  WARNING: directory {lib_key} not found, skipping.")
            continue

        for category in CATEGORIES:
            for panel_cfg in PANELS.values():
                for cell_line in panel_cfg["cell_lines"]:
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
                    df["Term_clean"] = (
                        df["Term"].str.replace(r"\s+", " ", regex=True).str.strip()
                    )
                    all_dfs.append(df)

    if not all_dfs:
        raise FileNotFoundError(
            f"No disease enrichment TSV files were found in {RESULTS_DIR}"
        )

    combined = pd.concat(all_dfs, ignore_index=True)
    print(f"Loaded {len(combined)} total disease enrichment entries")
    print(f"  Libraries: {combined['Library'].unique()}")
    print(f"  Categories: {combined['Category'].unique()}")
    print(f"  Cell lines: {combined['CellLine'].unique()}")
    return combined


def select_top_terms_by_library(df):
    """Pick top terms from each disease library for the current panel."""
    all_top_terms = []
    library_groups = []

    for lib_label in DISEASE_LIBRARIES.values():
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

        top_terms = list(dict.fromkeys(black_top + white_top))[:TOP_N_PER_LIB]
        if not top_terms:
            continue

        all_top_terms.extend(top_terms)
        library_groups.append((lib_label, top_terms))

    return all_top_terms, library_groups


def build_panel_tables(df, cell_lines, library_groups):
    """Build black and white matrices in the grouped row order."""
    all_terms = [term for _, terms in library_groups for term in terms]
    black_rows = []
    white_rows = []

    for lib_label, terms in library_groups:
        lib_df = df[df["Library"] == lib_label]
        for term in terms:
            for category, rows_list in [("black", black_rows), ("white", white_rows)]:
                cat_df = lib_df[
                    (lib_df["Category"] == category) & (lib_df["Term_clean"] == term)
                ]
                row = {}
                for cell_line in cell_lines:
                    val = cat_df.loc[cat_df["CellLine"] == cell_line, "neg_log10_adjp"]
                    row[cell_line] = val.values[0] if len(val) > 0 else 0.0
                rows_list.append(row)

    black_pivot = pd.DataFrame(black_rows, index=all_terms, columns=cell_lines)
    white_pivot = pd.DataFrame(white_rows, index=all_terms, columns=cell_lines)
    return black_pivot, white_pivot, all_terms


def plot_panel(df, panel_cfg):
    """Create one grouped disease heatmap panel."""
    cell_lines = panel_cfg["cell_lines"]
    panel_df = df[df["CellLine"].isin(cell_lines)].copy()

    all_terms, library_groups = select_top_terms_by_library(panel_df)
    if not all_terms:
        print(f"  WARNING: no terms available for panel {panel_cfg['subtitle']}, skipping.")
        return

    print(f"  Total terms across libraries: {len(all_terms)}")
    black_pivot, white_pivot, _ = build_panel_tables(panel_df, cell_lines, library_groups)

    display_columns = [panel_cfg["display_labels"].get(cl, cl) for cl in cell_lines]
    black_pivot.columns = display_columns
    white_pivot.columns = display_columns

    vmax = max(black_pivot.max().max(), white_pivot.max().max())
    vmax = max(vmax, 2)

    width_ratios = [max(2.2, 0.9 * len(cell_lines)), 5.4, max(2.2, 0.9 * len(cell_lines)), 0.3]
    fig = plt.figure(figsize=(5.5, 3.1))
    gs = fig.add_gridspec(1, 4, width_ratios=width_ratios, wspace=0.005)
    ax1 = fig.add_subplot(gs[0, 0])
    ax_labels = fig.add_subplot(gs[0, 1], sharey=ax1)
    ax2 = fig.add_subplot(gs[0, 2], sharey=ax1)
    cax = fig.add_subplot(gs[0, 3])

    cmap = sns.color_palette("YlOrRd", as_cmap=True) if sns is not None else plt.get_cmap("YlOrRd")

    draw_heatmap(
        ax=ax1,
        data=black_pivot,
        cmap=cmap,
        vmin=0,
        vmax=vmax,
    )
    ax1.set_title("Black", fontsize=TITLE_FS, fontweight="bold", pad=4)
    ax1.set_ylabel("")
    ax1.set_xlabel("")
    ax1.tick_params(axis="x", labelsize=XTICK_FS, rotation=35)
    for label in ax1.get_xticklabels():
        label.set_horizontalalignment("right")

    draw_heatmap(
        ax=ax2,
        data=white_pivot,
        cmap=cmap,
        vmin=0,
        vmax=vmax,
        cax=cax,
        cbar_label="-log₁₀(Adj. P)",
    )
    ax2.set_title("White", fontsize=TITLE_FS, fontweight="bold", pad=4)
    ax2.set_ylabel("")
    ax2.set_xlabel("")
    ax2.tick_params(axis="x", labelsize=XTICK_FS, rotation=35)
    for label in ax2.get_xticklabels():
        label.set_horizontalalignment("right")

    heatmap_ylim = ax1.get_ylim()
    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(heatmap_ylim)

    lib_colors = ["#f0f4ff", "#fff4f0", "#f0fff4", "#fffbe8"]
    row_idx = 0
    for grp_i, (lib_label, terms) in enumerate(library_groups):
        bg_color = lib_colors[grp_i % len(lib_colors)]
        ax_labels.axhspan(row_idx, row_idx + len(terms), color=bg_color, alpha=0.5)
        ax1.axhspan(row_idx, row_idx + len(terms), color=bg_color, alpha=0.15)
        ax2.axhspan(row_idx, row_idx + len(terms), color=bg_color, alpha=0.15)

        for i, term in enumerate(terms):
            ax_labels.text(
                0.5,
                row_idx + i + 0.5,
                format_term_two_lines(term),
                ha="center",
                va="center",
                fontsize=TERM_FS,
                fontweight="medium",
                linespacing=1.0,
            )

        row_idx += len(terms)

    ax_labels.axis("off")
    ax2.set_ylim(ax1.get_ylim())

    cax.tick_params(labelsize=CBAR_TICK_FS)
    cax.yaxis.label.set_size(CBAR_LABEL_FS)

    legend_patches = []
    for grp_i, (lib_label, _) in enumerate(library_groups):
        bg_color = lib_colors[grp_i % len(lib_colors)]
        legend_patches.append(
            mpatches.Patch(
                facecolor=bg_color,
                edgecolor="#999999",
                label=lib_label,
                alpha=0.7,
            )
        )

    fig.legend(
        handles=legend_patches,
        loc="lower center",
        ncol=len(library_groups),
        fontsize=LEGEND_FS,
        framealpha=0.9,
        edgecolor="#cccccc",
        bbox_to_anchor=(0.50, -0.08),
    )

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
        f"{panel_cfg['title']}\n{panel_cfg['subtitle']}",
        fontsize=SUPTITLE_FS,
        fontweight="bold",
        y=0.99,
    )

    plt.subplots_adjust(left=0.03, right=0.93, top=0.83, bottom=0.28)

    outpath_png = os.path.join(OUTPUT_DIR, f"{panel_cfg['output_stub']}.png")
    outpath_pdf = os.path.join(OUTPUT_DIR, f"{panel_cfg['output_stub']}.pdf")
    fig.savefig(outpath_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(outpath_pdf, bbox_inches="tight", facecolor="white")
    print(f"\nCombined disease heatmap saved: {outpath_png}")
    print(f"Combined disease heatmap saved: {outpath_pdf}")
    plt.close(fig)


def main():
    print("=" * 60)
    print("Combined Disease Enrichment — Black vs White Comparison")
    print("=" * 60)

    df = load_disease_results()

    for panel_name, panel_cfg in PANELS.items():
        print(
            f"\n--- Generating {panel_name} panel "
            f"(top {TOP_N_PER_LIB} per library) ---"
        )
        plot_panel(df, panel_cfg)

    print("\n✅ Done!")


if __name__ == "__main__":
    main()
