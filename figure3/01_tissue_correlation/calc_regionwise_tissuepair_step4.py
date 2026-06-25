"""
calc_regionwise_tissuepair_corr.py   (v3  •  Nature-Genetics-ready)
https://chatgpt.com/share/68af3820-0f28-8010-9ee7-59ca4227eba6
------------------------------------------------------------------
1.  region × tissue-pair Pearson-r matrix   → CSV
2.  publication-quality scatter-grid PDFs   → --pdf_dir

Example
-------
python calc_regionwise_tissuepair_corr.py \
       --meth       methylation_matrix.csv \
       --map        GTEx_811_SampleMap.txt \
       --out        region_tissuepair_corr.csv \
       --pdf_dir    region_scatter_grids \
       --min_pairs  10                     # 10 is now the default
"""

from __future__ import annotations
import argparse, itertools, re, sys, warnings
from pathlib  import Path
from typing   import Dict

import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ─────────────────────────── data helpers ───────────────────────────
def load_inputs(meth_fp: Path, map_fp: Path):
    meth = pd.read_csv(meth_fp, index_col=0)
    smap = pd.read_csv(map_fp, sep="\t")
    smap = smap[smap["SeqCore_ID"].isin(meth.columns)].reset_index(drop=True)
    meth = meth.loc[:, smap["SeqCore_ID"]]
    if meth.empty:
        sys.exit("✗  No overlapping sample IDs!")
    return meth, smap

def pick_one_sample_per_donor(smap: pd.DataFrame) -> Dict[str, Dict[str, str]]:
    """Return {tissue: {donor → chosen_sample}}  (first encounter kept)."""
    out = {}
    for tissue, grp in smap.groupby("Tissue"):
        out[tissue] = (
            grp.drop_duplicates("GTEx_ID")
               .set_index("GTEx_ID")["SeqCore_ID"]
               .to_dict()
        )
    return out

# ───────────────────────── correlations (unchanged) ─────────────────────────
def regionwise_corr(meth: pd.DataFrame,
                    by_tissue: Dict[str, Dict[str, str]],
                    min_pairs: int = 10) -> pd.DataFrame:
    tissues      = sorted(by_tissue)
    tissue_pairs = list(itertools.combinations(tissues, 2))
    out          = pd.DataFrame(index=meth.index, dtype=np.float32)

    for tA, tB in tissue_pairs:
        donors = list(set(by_tissue[tA]) & set(by_tissue[tB]))
        col    = np.full(len(meth), np.nan, dtype=np.float32)
        if len(donors) >= min_pairs:
            sA = meth[[by_tissue[tA][d] for d in donors]].to_numpy()
            sB = meth[[by_tissue[tB][d] for d in donors]].to_numpy()
            mask = (~np.isnan(sA)) & (~np.isnan(sB))
            enough = np.where(mask.sum(axis=1) >= min_pairs)[0]
            for i in enough:
                ai, bi = sA[i, mask[i]], sB[i, mask[i]]
                col[i] = np.corrcoef(ai, bi)[0, 1]
        out[f"{tA}‒{tB}"] = col
    return out

# ───────────────────────── publication scatter-grid ─────────────────────────
def safe_name(name: str) -> str:
    return re.sub(r"[^\w\-\.]", "_", name)

""" def scatter_grid(region: str,
                 meth: pd.DataFrame,
                 by_tissue: Dict[str, Dict[str, str]],
                 out_path: Path,
                 min_pairs: int = 10):
    tissues = sorted(by_tissue)
    # donors × tissues table (one chosen sample per donor/tissue)
    data = {}
    for t in tissues:
        s = by_tissue[t]
        col = meth.loc[region, list(s.values())].rename(
                    index={v: d for d, v in s.items()})
        data[t] = col
    df = pd.concat(data, axis=1).sort_index()

    if df.count().max() < min_pairs:          # not enough donors to show anything
        return

    # --- PairGrid with only lower-triangle scatter+regression lines ---
    sns.set_context("paper", font_scale=0.9)
    g = sns.PairGrid(df, vars=tissues, corner=True, height=2.2, aspect=1)
    for i in range(len(tissues)):
        ax = g.axes[i][i]          # main‐diagonal axis
        ax.set_visible(False)      # removes frame, ticks, labels

    def _reg_scatter(x, y, **kw):
        mask = ~np.isnan(x) & ~np.isnan(y)
        if mask.sum() < min_pairs:
            plt.gca().set_visible(False)
            return
        sns.regplot(x=x[mask], y=y[mask],
                    scatter_kws=dict(s=15, alpha=.8, edgecolor="none"),
                    line_kws=dict(lw=0.8), ci=None, truncate=False)
        # R² annotation
        r = np.corrcoef(x[mask], y[mask])[0, 1]
        plt.text(0.05, 0.9, f"$R^2$ = {r**2:.2f}",
                 transform=plt.gca().transAxes, fontsize=6)
        plt.xlim(0, 100)
        plt.ylim(0, 100)

    g.map_lower(_reg_scatter)

    # cosmetic: unify axis limits on all panels
    for ax in g.axes.flatten():
        if ax:
            ax.set_xlim(0, 100)
            ax.set_ylim(0, 100)

    g.fig.suptitle(region, y=1.015, fontsize=10)
    g.fig.tight_layout()
    g.fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(g.fig) """


def scatter_grid(region: str,
                 meth: pd.DataFrame,
                 by_tissue: Dict[str, Dict[str, str]],
                 out_path: Path,
                 min_pairs: int = 10):
    """
    Create a 1×5 grid of scatter plots: Y = each tissue, X = Whole Blood.
    Panels with < min_pairs donors are hidden. Axes fixed to 0–100. Adds R² & n.
    """
    import math
    sns.set_context("paper", font_scale=0.9)

    wb_name = "Whole Blood"
    tissues = sorted(by_tissue)
    if wb_name not in tissues:
        return

    # choose up to 5 non-WB tissues for a 1x5 layout
    others = [t for t in tissues if t != wb_name][:5]
    ncols = 5
    fig_w, fig_h = 2.3 * ncols, 2.3  # tuned for Nature Genetics aesthetics
    fig, axes = plt.subplots(1, ncols, figsize=(fig_w, fig_h))
    if ncols == 1:  # safety if ever ncols=1
        axes = [axes]

    for idx in range(ncols):
        ax = axes[idx]
        if idx >= len(others):
            ax.set_visible(False)
            continue

        t = others[idx]
        # donors with both tissues; align samples donor-wise
        donors = sorted(set(by_tissue[t]).intersection(by_tissue[wb_name]))
        if not donors:
            ax.set_visible(False)
            continue

        x_ids = [by_tissue[wb_name][d] for d in donors]  # Whole Blood
        y_ids = [by_tissue[t][d]        for d in donors]  # Tissue t

        x = meth.loc[region, x_ids].to_numpy(dtype=float)
        y = meth.loc[region, y_ids].to_numpy(dtype=float)

        mask = (~np.isnan(x)) & (~np.isnan(y))
        n_ok = int(mask.sum())
        if n_ok < min_pairs:
            ax.set_visible(False)
            continue

        # scatter + regression; fixed axes; small, sharp markers
        sns.regplot(x=x[mask], y=y[mask], ax=ax,
                    scatter_kws=dict(s=15, alpha=.85, edgecolor="none"),
                    line_kws=dict(lw=1.0), ci=None, truncate=False)
        r = np.corrcoef(x[mask], y[mask])[0, 1]
        ax.text(0.04, 0.92, f"$R^2$ = {r**2:.2f}\n$n$ = {n_ok}",
                transform=ax.transAxes, fontsize=6)

        ax.set_xlim(0, 100); ax.set_ylim(0, 100)
        ax.set_xlabel(f"{wb_name} (%)", fontsize=8)
        ax.set_ylabel(f"{t} (%)", fontsize=8)
        ax.tick_params(labelsize=7)

    fig.suptitle(region, y=0.995, fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ─────────────────────────────── main ──────────────────────────────
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--meth",       required=True, type=Path)
    p.add_argument("--map",        required=True, type=Path)
    p.add_argument("--out",        required=True, type=Path,
                   help="CSV filename for correlation matrix")
    p.add_argument("--pdf_dir",    type=Path,
                   help="Directory for one PDF per region (optional)")
    p.add_argument("--min_pairs",  default=10, type=int,
                   help="Min donor pairs required for r/R² (default 10)")
    p.add_argument("--max_regions", default=None, type=int,
                   help="Limit #regions processed for PDFs (speed test)")
    args = p.parse_args()

    sns.set_style("white")  # clean, publication-style
    warnings.filterwarnings("ignore", category=UserWarning)

    meth, smap  = load_inputs(args.meth, args.map)
    by_tissue   = pick_one_sample_per_donor(smap)

    # 1.  correlations → CSV
    corr = regionwise_corr(meth, by_tissue, args.min_pairs)
    corr.to_csv(args.out)
    print(f"✓  correlation matrix  {args.out}  "
          f"({corr.shape[0]:,} regions × {corr.shape[1]} pairs)")

    # 2.  scatter-grid PDFs
    if args.pdf_dir:
        outdir = args.pdf_dir
        outdir.mkdir(parents=True, exist_ok=True)
        regions = meth.index if args.max_regions is None else meth.index[:args.max_regions]
        for region in regions:
            pdf = outdir / f"{safe_name(region)}.pdf"
            scatter_grid(region, meth, by_tissue, pdf, args.min_pairs)
        n = len(list(outdir.glob("*.pdf")))
        print(f"✓  scatter-grid PDFs  →  {outdir}  ({n} files)")

if __name__ == "__main__":
    main()

