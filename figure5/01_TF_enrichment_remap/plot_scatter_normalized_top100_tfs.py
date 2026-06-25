import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
from adjustText import adjust_text

# ── Column names (BED file basenames used as query_file values) ───────────────
WHITE_COL = "White_CpGge5_Rangegt20.bed"
BLACK_COL = "Black_CpGge5_Rangegt20.bed"
RAND_COL  = "random_background.bed"

# Cell-line display order — fills 2×3 grid row by row; slot [1,2] stays empty
CELL_LINE_ORDER = ["HepG2", "K562", "MCF7", "WA01", "WA09"]

# ── Font sizes ────────────────────────────────────────────────────────────────
FS_TITLE   = 20   # panel title (cell line name)
FS_AXLABEL = 18   # axis labels
FS_TICK    = 16   # tick labels
FS_TFLABEL = 14   # TF text labels (many per panel → keep small but legible)
FS_LEGEND  = 16   # legend

# ── Grid layout ───────────────────────────────────────────────────────────────
N_ROWS, N_COLS = 2, 3

# ── Output ────────────────────────────────────────────────────────────────────
outdir = Path("overlap_results_norm/scatter_plots_norm")
outdir.mkdir(parents=True, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
df = pd.read_csv("overlap_results_norm/tf_overlap_normalized_counts.tsv", sep="\t")

# ── Build per-cell-line matrices (top 100 TFs, zero-filled) ──────────────────
mats = {}
for cell_line in CELL_LINE_ORDER:
    sub = df[df["cell_line"] == cell_line]
    if sub.empty:
        print(f"[WARNING] No data for {cell_line}, panel will be blank.")
        mats[cell_line] = None
        continue

    # Pivot: TF × query_file; fillna(0) → TFs absent from a query column → 0
    mat = (
        sub.groupby(["tf", "query_file"], as_index=False)["normalized_count"]
           .sum()
           .pivot(index="tf", columns="query_file", values="normalized_count")
           .fillna(0)   # TFs missing from a query series get 0, not NaN
    )

    # Ensure all three query columns exist (edge case: no overlaps at all)
    for col in [WHITE_COL, BLACK_COL, RAND_COL]:
        if col not in mat.columns:
            mat[col] = 0.0

    # Diagnostics
    n_bz = (mat[BLACK_COL] == 0).sum()
    n_wz = (mat[WHITE_COL] == 0).sum()
    print(f"[{cell_line}] {len(mat)} TFs total  — Black=0: {n_bz}, White=0: {n_wz}")

    # Top 100 by total normalized count; zero-filled TFs are still eligible
    top100 = mat.sum(axis=1).sort_values(ascending=False).head(100).index
    mat    = mat.loc[top100]

    n_bz2 = (mat[BLACK_COL] == 0).sum()
    n_wz2 = (mat[WHITE_COL] == 0).sum()
    print(f"  After top-100 cut → Black=0: {n_bz2}, White=0: {n_wz2}")

    mats[cell_line] = mat

# ── Create figure — portrait letter (8.5 × 11 in), 2 rows × 3 cols ───────────
# Each panel uses equal aspect so the axes box itself is square.
# Since x and y are both "norm. overlaps / 1k" the equal aspect is also
# scientifically meaningful (the y=x diagonal becomes a true 45° line).
fig = plt.figure(figsize=(8.5, 9.0))
gs  = gridspec.GridSpec(
    N_ROWS, N_COLS,
    figure=fig,
    left=0.09, right=0.97,
    top=0.93, bottom=0.07,
    wspace=0.40,   # horizontal gap between panels
    hspace=0.45,   # vertical gap between rows
)

for panel_idx, cell_line in enumerate(CELL_LINE_ORDER):
    row = panel_idx // N_COLS
    col = panel_idx %  N_COLS
    ax  = fig.add_subplot(gs[row, col])
    mat = mats.get(cell_line)

    if mat is None:
        ax.set_visible(False)
        continue

    # ── Scatter points ────────────────────────────────────────────────────────
    ax.scatter(
        mat[WHITE_COL], mat[BLACK_COL],
        color="black", alpha=0.75, s=25, zorder=3,
        label="Black CoRSIVs"
    )
    ax.scatter(
        mat[WHITE_COL], mat[RAND_COL],
        color="hotpink", alpha=0.75, s=25, zorder=3,
        label="Random"
    )

    # ── y = x reference diagonal ──────────────────────────────────────────────
    max_val = max(mat[WHITE_COL].max(), mat[BLACK_COL].max(), mat[RAND_COL].max())
    if max_val > 0:
        ax.plot([0, max_val], [0, max_val], "k--", lw=0.9, alpha=0.35, label="y = x")

    # ── TF text labels ────────────────────────────────────────────────────────
    texts = []
    top20_tfs = mat[WHITE_COL].sort_values(ascending=False).head(20).index
    for tf in top20_tfs:
        x  = mat.loc[tf, WHITE_COL]
        yb = mat.loc[tf, BLACK_COL]
        yr = mat.loc[tf, RAND_COL]
        texts.append(
            ax.text(x, yb, tf, fontsize=FS_TFLABEL, color="black",
                    ha="center", va="center")
        )

    # ── adjustText: repel labels from each other and from data points ─────────
    adjust_text(
        texts,
        ax=ax,
        expand=(1.2, 1.4),
        force_text=(0.10, 0.15),
        force_points=(0.15, 0.20),
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.4, alpha=0.55),
        time_lim=8,              # allow up to 8 s per panel for convergence
    )

    # ── Make axes box square (equal data units on x and y) ───────────────────
    ax.set_aspect("equal", adjustable="box")

    # ── Axes labels & formatting ──────────────────────────────────────────────
    ax.set_title(cell_line, fontsize=FS_TITLE, fontweight="bold", pad=7)
    ax.set_xlabel("White (TF-bound regions / 1k regions)", fontsize=FS_AXLABEL, labelpad=4)
    # Show y-label only on left column panels
    if col == 0:
        ax.set_ylabel("Black / Random\n(TF-bound regions / 1k regions)",
                      fontsize=FS_AXLABEL, labelpad=4)
    else:
        ax.set_ylabel("")

    ax.tick_params(axis="both", labelsize=FS_TICK)

    # Legend in the top-left corner of the first panel only
    if panel_idx == 0:
        ax.legend(fontsize=FS_LEGEND, framealpha=0.75,
                  loc="upper left", handletextpad=0.4)

# Hide the unused 6th subplot slot
unused_row = (len(CELL_LINE_ORDER)) // N_COLS
unused_col = (len(CELL_LINE_ORDER)) %  N_COLS
fig.add_subplot(gs[unused_row, unused_col]).set_visible(False)

# ── Super-title ───────────────────────────────────────────────────────────────
fig.suptitle(
    "TF Overlaps per 1,000 Regions — Top 100 TFs  (Black CoRSIVs & Random vs White)",
    fontsize=24, fontweight="bold", y=0.97
)

# ── Save grid figure ──────────────────────────────────────────────────────────
pdf_path = outdir / "scatter_norm_top100_all_panels.pdf"
fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
print(f"\nGrid PDF saved : {pdf_path.absolute()}")

png_path = outdir / "scatter_norm_top100_all_panels.png"
fig.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
print(f"Grid PNG saved : {png_path.absolute()}")

plt.close(fig)

# ══════════════════════════════════════════════════════════════════════════════
# Individual per-cell-line publication PDFs
# Same style as the grid panels but on a larger canvas so fonts and labels
# can breathe. Each output is a square 7×7-inch PDF.
# ══════════════════════════════════════════════════════════════════════════════

# Font sizes for the larger individual canvas
IFS_TITLE   = 24
IFS_AXLABEL = 22
IFS_TICK    = 20
IFS_TFLABEL = 18
IFS_LEGEND  = 20

print("\n── Individual cell-line PDFs ─────────────────────────────────────────────")
for cell_line in CELL_LINE_ORDER:
    mat = mats.get(cell_line)
    if mat is None:
        print(f"  [SKIP] {cell_line} — no data")
        continue

    fig_s, ax_s = plt.subplots(figsize=(7, 7))

    # ── Scatter points ────────────────────────────────────────────────────────
    ax_s.scatter(
        mat[WHITE_COL], mat[BLACK_COL],
        color="black", alpha=0.75, s=35, zorder=3,
        label="Black CoRSIVs"
    )
    ax_s.scatter(
        mat[WHITE_COL], mat[RAND_COL],
        color="hotpink", alpha=0.75, s=35, zorder=3,
        label="Random"
    )

    # ── y = x reference diagonal ──────────────────────────────────────────────
    max_val = max(mat[WHITE_COL].max(), mat[BLACK_COL].max(), mat[RAND_COL].max())
    if max_val > 0:
        ax_s.plot([0, max_val], [0, max_val], "k--", lw=1.0, alpha=0.35, label="y = x")

    # ── TF text labels ────────────────────────────────────────────────────────
    texts_s = []
    top20_tfs_s = mat[WHITE_COL].sort_values(ascending=False).head(20).index
    for tf in top20_tfs_s:
        x  = mat.loc[tf, WHITE_COL]
        yb = mat.loc[tf, BLACK_COL]
        yr = mat.loc[tf, RAND_COL]
        texts_s.append(
            ax_s.text(x, yb, tf, fontsize=IFS_TFLABEL, color="black",
                      ha="center", va="center")
        )

    # ── adjustText repulsion ──────────────────────────────────────────────────
    adjust_text(
        texts_s,
        ax=ax_s,
        expand=(1.2, 1.4),
        force_text=(0.10, 0.15),
        force_points=(0.15, 0.20),
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.5, alpha=0.55),
        time_lim=12,       # more time on the larger canvas for better convergence
    )

    # ── Square axes (equal data units on both axes) ───────────────────────────
    ax_s.set_aspect("equal", adjustable="box")

    # ── Labels & formatting ───────────────────────────────────────────────────
    ax_s.set_title(
        f"{cell_line}  —  TF Overlaps per 1,000 Regions (Top 100 TFs)",
        fontsize=IFS_TITLE, fontweight="bold", pad=10
    )
    ax_s.set_xlabel("White (TF-bound regions / 1k regions)", fontsize=IFS_AXLABEL, labelpad=6)
    ax_s.set_ylabel("Black / Random (TF-bound regions / 1k regions)", fontsize=IFS_AXLABEL, labelpad=6)
    ax_s.tick_params(axis="both", labelsize=IFS_TICK)
    ax_s.legend(fontsize=IFS_LEGEND, framealpha=0.75,
                loc="upper left", handletextpad=0.4)

    fig_s.tight_layout()

    # ── Save individual PDF + PNG ─────────────────────────────────────────────
    ind_pdf = outdir / f"scatter_norm_top100_{cell_line}.pdf"
    fig_s.savefig(ind_pdf, format="pdf", bbox_inches="tight")

    ind_png = outdir / f"scatter_norm_top100_{cell_line}.png"
    fig_s.savefig(ind_png, format="png", dpi=300, bbox_inches="tight")

    plt.close(fig_s)
    print(f"  {cell_line}: PDF + PNG saved → {ind_pdf.name}")
