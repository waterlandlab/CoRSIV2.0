import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
from adjustText import adjust_text

# BED basenames used in overlap_results_norm/tf_overlap_normalized_counts.tsv
WHITE_RANDOM_COL = "White_random_background.bed"
BLACK_RANDOM_COL = "Black_random_background.bed"

CELL_LINE_ORDER = ["HepG2", "K562", "MCF7", "WA01", "WA09"]

FS_TITLE = 20
FS_AXLABEL = 18
FS_TICK = 16
FS_TFLABEL = 14
FS_LEGEND = 16

N_ROWS, N_COLS = 2, 3

outdir = Path("overlap_results_norm/scatter_plots_random_sets")
outdir.mkdir(parents=True, exist_ok=True)

counts_path = Path("overlap_results_norm/tf_overlap_normalized_counts.tsv")
if not counts_path.exists():
    raise SystemExit(f"Missing normalized counts table: {counts_path}")

df = pd.read_csv(counts_path, sep="\t")
present_queries = set(df["query_file"].dropna().unique())
needed_queries = {WHITE_RANDOM_COL, BLACK_RANDOM_COL}
missing_queries = sorted(needed_queries - present_queries)
if missing_queries:
    missing_str = ", ".join(missing_queries)
    raise SystemExit(
        "The normalized counts table does not include the random-set query files: "
        f"{missing_str}\n"
        "Re-run count_overlap_normalized_tfs.sh to regenerate "
        "overlap_results_norm/tf_overlap_normalized_counts.tsv."
    )

mats = {}
for cell_line in CELL_LINE_ORDER:
    sub = df[
        (df["cell_line"] == cell_line)
        & (df["query_file"].isin([WHITE_RANDOM_COL, BLACK_RANDOM_COL]))
    ]
    if sub.empty:
        print(f"[WARNING] No data for {cell_line}, panel will be blank.")
        mats[cell_line] = None
        continue

    mat = (
        sub.groupby(["tf", "query_file"], as_index=False)["normalized_count"]
        .sum()
        .pivot(index="tf", columns="query_file", values="normalized_count")
        .fillna(0)
    )

    for col in [WHITE_RANDOM_COL, BLACK_RANDOM_COL]:
        if col not in mat.columns:
            mat[col] = 0.0

    n_white_zero = (mat[WHITE_RANDOM_COL] == 0).sum()
    n_black_zero = (mat[BLACK_RANDOM_COL] == 0).sum()
    print(
        f"[{cell_line}] {len(mat)} TFs total  — "
        f"White-random=0: {n_white_zero}, Black-random=0: {n_black_zero}"
    )

    top100 = mat.sum(axis=1).sort_values(ascending=False).head(100).index
    mat = mat.loc[top100]

    n_white_zero_top = (mat[WHITE_RANDOM_COL] == 0).sum()
    n_black_zero_top = (mat[BLACK_RANDOM_COL] == 0).sum()
    print(
        f"  After top-100 cut → White-random=0: {n_white_zero_top}, "
        f"Black-random=0: {n_black_zero_top}"
    )

    mats[cell_line] = mat

fig = plt.figure(figsize=(8.5, 9.0))
gs = gridspec.GridSpec(
    N_ROWS,
    N_COLS,
    figure=fig,
    left=0.09,
    right=0.97,
    top=0.93,
    bottom=0.07,
    wspace=0.40,
    hspace=0.45,
)

for panel_idx, cell_line in enumerate(CELL_LINE_ORDER):
    row = panel_idx // N_COLS
    col = panel_idx % N_COLS
    ax = fig.add_subplot(gs[row, col])
    mat = mats.get(cell_line)

    if mat is None:
        ax.set_visible(False)
        continue

    ax.scatter(
        mat[WHITE_RANDOM_COL],
        mat[BLACK_RANDOM_COL],
        color="steelblue",
        alpha=0.80,
        s=28,
        zorder=3,
        label="TF overlap"
    )

    max_val = max(mat[WHITE_RANDOM_COL].max(), mat[BLACK_RANDOM_COL].max())
    if max_val > 0:
        ax.plot([0, max_val], [0, max_val], "k--", lw=0.9, alpha=0.35, label="y = x")

    texts = []
    top20_tfs = mat[WHITE_RANDOM_COL].sort_values(ascending=False).head(20).index
    for tf in top20_tfs:
        x = mat.loc[tf, WHITE_RANDOM_COL]
        y = mat.loc[tf, BLACK_RANDOM_COL]
        texts.append(
            ax.text(
                x,
                y,
                tf,
                fontsize=FS_TFLABEL,
                color="black",
                ha="center",
                va="center",
            )
        )

    adjust_text(
        texts,
        ax=ax,
        expand=(1.2, 1.4),
        force_text=(0.10, 0.15),
        force_points=(0.15, 0.20),
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.4, alpha=0.55),
        time_lim=8,
    )

    ax.set_aspect("equal", adjustable="box")
    ax.set_title(cell_line, fontsize=FS_TITLE, fontweight="bold", pad=7)
    ax.set_xlabel(
        "White random\n(TF-bound regions / 1k regions)",
        fontsize=FS_AXLABEL,
        labelpad=4,
    )
    if col == 0:
        ax.set_ylabel(
            "Black random\n(TF-bound regions / 1k regions)",
            fontsize=FS_AXLABEL,
            labelpad=4,
        )
    else:
        ax.set_ylabel("")
    ax.tick_params(axis="both", labelsize=FS_TICK)

    if panel_idx == 0:
        ax.legend(fontsize=FS_LEGEND, framealpha=0.75, loc="upper left", handletextpad=0.4)

unused_row = len(CELL_LINE_ORDER) // N_COLS
unused_col = len(CELL_LINE_ORDER) % N_COLS
fig.add_subplot(gs[unused_row, unused_col]).set_visible(False)

fig.suptitle(
    "TF Overlaps per 1,000 Regions — Top 100 TFs  (Black-random vs White-random)",
    fontsize=24,
    fontweight="bold",
    y=0.97,
)

grid_pdf = outdir / "scatter_norm_top100_random_sets_all_panels.pdf"
fig.savefig(grid_pdf, format="pdf", bbox_inches="tight")
print(f"\nGrid PDF saved : {grid_pdf.absolute()}")

grid_png = outdir / "scatter_norm_top100_random_sets_all_panels.png"
fig.savefig(grid_png, format="png", dpi=300, bbox_inches="tight")
print(f"Grid PNG saved : {grid_png.absolute()}")

plt.close(fig)

IFS_TITLE = 24
IFS_AXLABEL = 22
IFS_TICK = 20
IFS_TFLABEL = 18
IFS_LEGEND = 20

print("\n-- Individual cell-line PDFs --------------------------------------------")
for cell_line in CELL_LINE_ORDER:
    mat = mats.get(cell_line)
    if mat is None:
        print(f"  [SKIP] {cell_line} — no data")
        continue

    fig_s, ax_s = plt.subplots(figsize=(7, 7))

    ax_s.scatter(
        mat[WHITE_RANDOM_COL],
        mat[BLACK_RANDOM_COL],
        color="steelblue",
        alpha=0.80,
        s=38,
        zorder=3,
        label="TF overlap"
    )

    max_val = max(mat[WHITE_RANDOM_COL].max(), mat[BLACK_RANDOM_COL].max())
    if max_val > 0:
        ax_s.plot([0, max_val], [0, max_val], "k--", lw=1.0, alpha=0.35, label="y = x")

    texts_s = []
    top20_tfs_s = mat[WHITE_RANDOM_COL].sort_values(ascending=False).head(20).index
    for tf in top20_tfs_s:
        x = mat.loc[tf, WHITE_RANDOM_COL]
        y = mat.loc[tf, BLACK_RANDOM_COL]
        texts_s.append(
            ax_s.text(
                x,
                y,
                tf,
                fontsize=IFS_TFLABEL,
                color="black",
                ha="center",
                va="center",
            )
        )

    adjust_text(
        texts_s,
        ax=ax_s,
        expand=(1.2, 1.4),
        force_text=(0.10, 0.15),
        force_points=(0.15, 0.20),
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.5, alpha=0.55),
        time_lim=12,
    )

    ax_s.set_aspect("equal", adjustable="box")
    ax_s.set_title(
        f"{cell_line}  —  TF Overlaps per 1,000 Regions (Top 100 TFs)",
        fontsize=IFS_TITLE,
        fontweight="bold",
        pad=10,
    )
    ax_s.set_xlabel(
        "White random (TF-bound regions / 1k regions)",
        fontsize=IFS_AXLABEL,
        labelpad=6,
    )
    ax_s.set_ylabel(
        "Black random (TF-bound regions / 1k regions)",
        fontsize=IFS_AXLABEL,
        labelpad=6,
    )
    ax_s.tick_params(axis="both", labelsize=IFS_TICK)
    ax_s.legend(fontsize=IFS_LEGEND, framealpha=0.75, loc="upper left", handletextpad=0.4)

    fig_s.tight_layout()

    ind_pdf = outdir / f"scatter_norm_top100_random_sets_{cell_line}.pdf"
    fig_s.savefig(ind_pdf, format="pdf", bbox_inches="tight")

    ind_png = outdir / f"scatter_norm_top100_random_sets_{cell_line}.png"
    fig_s.savefig(ind_png, format="png", dpi=300, bbox_inches="tight")

    plt.close(fig_s)
