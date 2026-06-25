import io
import os
import subprocess
from collections import Counter
from contextlib import redirect_stdout
from pathlib import Path

import pandas as pd

OUTDIR = Path("overlap_results_norm/scatter_plots_matched_controls")
OUTDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUTDIR / ".mplconfig"))
(OUTDIR / ".mplconfig").mkdir(parents=True, exist_ok=True)

WHITE_MATCHED_COL = "white.sorted.matched.controls.bed"
BLACK_MATCHED_COL = "black.sorted.matched.controls.bed"

CELL_LINE_ORDER = ["HepG2", "K562", "MCF7", "WA01", "WA09"]

FS_TITLE = 20
FS_AXLABEL = 18
FS_TICK = 16
FS_TFLABEL = 14
FS_LEGEND = 16

IFS_TITLE = 24
IFS_AXLABEL = 22
IFS_TICK = 20
IFS_TFLABEL = 18
IFS_LEGEND = 20

N_ROWS, N_COLS = 2, 3

query_dir = Path("query")
remap_glob = "remap2022_all_macs2_hg38_v1_0_*.bed"
outdir = OUTDIR

from adjustText import adjust_text
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

counts_path = outdir / "tf_overlap_normalized_counts_matched_controls.tsv"


def count_nonempty_lines(path: Path) -> int:
    with path.open() as handle:
        return sum(1 for line in handle if line.strip())


def count_query_columns(path: Path) -> int:
    with path.open() as handle:
        for line in handle:
            if line.strip():
                return len(line.rstrip("\n").split("\t"))
    raise ValueError(f"Query BED is empty: {path}")


def extract_tf(remap_name_field: str) -> str:
    parts = remap_name_field.split(".")
    return parts[1] if len(parts) >= 2 else remap_name_field


def compute_tf_counts(query_path: Path, remap_path: Path, query_ncols: int) -> Counter:
    cmd = [
        "bedtools",
        "intersect",
        "-a",
        str(query_path),
        "-b",
        str(remap_path),
        "-wa",
        "-wb",
    ]
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    counts = Counter()
    seen = set()

    assert proc.stdout is not None
    for line in proc.stdout:
        cols = line.rstrip("\n").split("\t")
        if len(cols) <= query_ncols + 3:
            continue

        tf = extract_tf(cols[query_ncols + 3])
        qkey = ":".join(cols[:3])
        seen_key = (tf, qkey)
        if seen_key in seen:
            continue
        seen.add(seen_key)
        counts[tf] += 1

    stderr = ""
    if proc.stderr is not None:
        stderr = proc.stderr.read().strip()
    return_code = proc.wait()
    if return_code != 0:
        raise RuntimeError(
            f"bedtools intersect failed for {query_path.name} vs {remap_path.name}: {stderr}"
        )

    return counts


def build_counts_table() -> pd.DataFrame:
    query_files = [
        query_dir / WHITE_MATCHED_COL,
        query_dir / BLACK_MATCHED_COL,
    ]
    remap_files = sorted(Path(".").glob(remap_glob))

    missing_queries = [str(path) for path in query_files if not path.exists()]
    if missing_queries:
        raise SystemExit(
            "Missing matched-control query BED file(s): " + ", ".join(missing_queries)
        )
    if not remap_files:
        raise SystemExit("No remap BED files found in the current directory.")

    rows = []
    for query_path in query_files:
        total_regions = count_nonempty_lines(query_path)
        if total_regions == 0:
            raise SystemExit(f"Query BED is empty: {query_path}")

        query_ncols = count_query_columns(query_path)
        print(f"[INFO] Processing {query_path.name} ({total_regions} regions)")

        for remap_path in remap_files:
            cell_line = remap_path.stem.replace("remap2022_all_macs2_hg38_v1_0_", "")
            tf_counts = compute_tf_counts(query_path, remap_path, query_ncols)
            print(
                f"  {remap_path.name}: {len(tf_counts)} TFs with overlaps in {cell_line}"
            )

            for tf, query_count in sorted(tf_counts.items()):
                rows.append(
                    {
                        "query_file": query_path.name,
                        "remap_file": remap_path.name,
                        "cell_line": cell_line,
                        "tf": tf,
                        "query_count": query_count,
                        "normalized_count": (query_count / total_regions) * 1000,
                    }
                )

    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("No overlaps detected for the matched-control query files.")

    df.to_csv(counts_path, sep="\t", index=False)
    print(f"\nMatched-control counts written: {counts_path.absolute()}")
    return df


def build_matrices(df: pd.DataFrame) -> dict:
    mats = {}
    for cell_line in CELL_LINE_ORDER:
        sub = df[
            (df["cell_line"] == cell_line)
            & (df["query_file"].isin([WHITE_MATCHED_COL, BLACK_MATCHED_COL]))
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

        for col in [WHITE_MATCHED_COL, BLACK_MATCHED_COL]:
            if col not in mat.columns:
                mat[col] = 0.0

        n_white_zero = (mat[WHITE_MATCHED_COL] == 0).sum()
        n_black_zero = (mat[BLACK_MATCHED_COL] == 0).sum()
        print(
            f"[{cell_line}] {len(mat)} TFs total  - "
            f"White-matched=0: {n_white_zero}, Black-matched=0: {n_black_zero}"
        )

        top100 = mat.sum(axis=1).sort_values(ascending=False).head(100).index
        mat = mat.loc[top100]

        n_white_zero_top = (mat[WHITE_MATCHED_COL] == 0).sum()
        n_black_zero_top = (mat[BLACK_MATCHED_COL] == 0).sum()
        print(
            f"  After top-100 cut -> White-matched=0: {n_white_zero_top}, "
            f"Black-matched=0: {n_black_zero_top}"
        )

        mats[cell_line] = mat

    return mats


def add_scatter(ax, mat: pd.DataFrame, point_size: int, tf_fontsize: int, line_width: float, time_limit: int):
    ax.scatter(
        mat[WHITE_MATCHED_COL],
        mat[BLACK_MATCHED_COL],
        color="steelblue",
        alpha=0.80,
        s=point_size,
        zorder=3,
        label="TF overlap",
    )

    max_val = max(mat[WHITE_MATCHED_COL].max(), mat[BLACK_MATCHED_COL].max())
    if max_val > 0:
        ax.plot([0, max_val], [0, max_val], "k--", lw=line_width, alpha=0.35, label="y = x")

    texts = []
    top20_tfs = mat[WHITE_MATCHED_COL].sort_values(ascending=False).head(20).index
    for tf in top20_tfs:
        x = mat.loc[tf, WHITE_MATCHED_COL]
        y = mat.loc[tf, BLACK_MATCHED_COL]
        texts.append(
            ax.text(
                x,
                y,
                tf,
                fontsize=tf_fontsize,
                color="black",
                ha="center",
                va="center",
            )
        )

    with redirect_stdout(io.StringIO()):
        adjust_text(
            texts,
            ax=ax,
            expand=(1.2, 1.4),
            force_text=(0.10, 0.15),
            force_points=(0.15, 0.20),
            arrowprops=dict(
                arrowstyle="-",
                color="gray",
                lw=0.4 if line_width < 1 else 0.5,
                alpha=0.55,
            ),
            time_lim=time_limit,
        )


def plot_grid(mats: dict) -> None:
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

        add_scatter(ax, mat, point_size=28, tf_fontsize=FS_TFLABEL, line_width=0.9, time_limit=8)

        ax.set_aspect("equal", adjustable="box")
        ax.set_title(cell_line, fontsize=FS_TITLE, fontweight="bold", pad=7)
        ax.set_xlabel(
            "White matched controls\n(TF-bound regions / 1k regions)",
            fontsize=FS_AXLABEL,
            labelpad=4,
        )
        if col == 0:
            ax.set_ylabel(
                "Black matched controls\n(TF-bound regions / 1k regions)",
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
        "TF Overlaps per 1,000 Regions - Top 100 TFs  (Black-matched vs White-matched)",
        fontsize=24,
        fontweight="bold",
        y=0.97,
    )

    grid_pdf = outdir / "scatter_norm_top100_matched_controls_all_panels.pdf"
    fig.savefig(grid_pdf, format="pdf", bbox_inches="tight")
    print(f"Grid PDF saved: {grid_pdf.absolute()}")

    grid_png = outdir / "scatter_norm_top100_matched_controls_all_panels.png"
    fig.savefig(grid_png, format="png", dpi=300, bbox_inches="tight")
    print(f"Grid PNG saved: {grid_png.absolute()}")

    plt.close(fig)


def plot_individual_panels(mats: dict) -> None:
    print("\n-- Individual cell-line plots -------------------------------------------")
    for cell_line in CELL_LINE_ORDER:
        mat = mats.get(cell_line)
        if mat is None:
            print(f"  [SKIP] {cell_line} - no data")
            continue

        fig_s, ax_s = plt.subplots(figsize=(7, 7))
        add_scatter(
            ax_s,
            mat,
            point_size=38,
            tf_fontsize=IFS_TFLABEL,
            line_width=1.0,
            time_limit=12,
        )

        ax_s.set_aspect("equal", adjustable="box")
        ax_s.set_title(
            f"{cell_line}  -  TF Overlaps per 1,000 Regions (Top 100 TFs)",
            fontsize=IFS_TITLE,
            fontweight="bold",
            pad=10,
        )
        ax_s.set_xlabel(
            "White matched controls (TF-bound regions / 1k regions)",
            fontsize=IFS_AXLABEL,
            labelpad=6,
        )
        ax_s.set_ylabel(
            "Black matched controls (TF-bound regions / 1k regions)",
            fontsize=IFS_AXLABEL,
            labelpad=6,
        )
        ax_s.tick_params(axis="both", labelsize=IFS_TICK)
        ax_s.legend(fontsize=IFS_LEGEND, framealpha=0.75, loc="upper left", handletextpad=0.4)

        pdf_path = outdir / f"scatter_norm_top100_matched_controls_{cell_line}.pdf"
        png_path = outdir / f"scatter_norm_top100_matched_controls_{cell_line}.png"
        fig_s.savefig(pdf_path, format="pdf", bbox_inches="tight")
        fig_s.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
        plt.close(fig_s)
        print(f"  Saved {cell_line}: {pdf_path.name}, {png_path.name}")


def main() -> None:
    df = build_counts_table()
    mats = build_matrices(df)
    plot_grid(mats)
    plot_individual_panels(mats)


if __name__ == "__main__":
    main()
