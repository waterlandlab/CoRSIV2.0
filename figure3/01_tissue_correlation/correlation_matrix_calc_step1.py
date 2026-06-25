"""
calc_regionwise_tissuepair_corr.py
----------------------------------
Run locally with:
    python calc_regionwise_tissuepair_corr.py \
        --meth  methylation_matrix.csv \
        --map   GTEx_811_SampleMap.txt \
        --min_pairs 3 \
        --out  region_tissuepair_corr.csv

Outputs a wide matrix:
    index  = region_id
    columns= "<TissueA>‒<TissueB>"
    cells  = Pearson r   (NaN if fewer than --min_pairs donors)

Author: C J’s ChatGPT helper
"""
from __future__ import annotations

import argparse
import itertools
import numpy  as np
import pandas as pd
from pathlib import Path
from typing  import Dict, List

# ----------------------------------------------------------------------
def load_inputs(meth_fp: Path, map_fp: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load methylation matrix and sample map; keep only intersecting samples."""
    meth = pd.read_csv(meth_fp, index_col=0)
    smap = pd.read_csv(map_fp, sep="\t")

    common = smap["SeqCore_ID"].isin(meth.columns)
    smap   = smap[common].copy()
    meth   = meth.loc[:, smap["SeqCore_ID"]]

    if meth.shape[1] == 0:
        raise RuntimeError("No sample IDs overlap between matrix and sample-map!")

    return meth, smap.reset_index(drop=True)

# ----------------------------------------------------------------------
def pick_one_sample_per_donor(smap: pd.DataFrame) -> Dict[str, Dict[str, str]]:
    """
    Build per-tissue dictionaries {donor → chosen_sample_id}.
    If a donor has >1 sample for a tissue we keep the FIRST occurrence
    (you could change this to `.mean()` or `.median()` on the methylation
    values—see comment in main()).
    """
    by_tissue: Dict[str, Dict[str, str]] = {}
    for tissue, grp in smap.groupby("Tissue"):
        # drop_duplicates keeps first row per GTEx_ID
        by_tissue[tissue] = (
            grp.drop_duplicates("GTEx_ID")
               .set_index("GTEx_ID")["SeqCore_ID"]
               .to_dict()
        )
    return by_tissue

# ----------------------------------------------------------------------
def regionwise_corr(meth: pd.DataFrame,
                    by_tissue: Dict[str, Dict[str, str]],
                    min_pairs: int = 3) -> pd.DataFrame:
    """
    Compute Pearson r for every region × tissue-pair, requiring ≥min_pairs donors.
    """
    tissues     = sorted(by_tissue.keys())
    tissue_pairs= list(itertools.combinations(tissues, 2))
    result      = pd.DataFrame(index=meth.index, dtype=np.float32)

    for tA, tB in tissue_pairs:
        donors = list(set(by_tissue[tA]).intersection(by_tissue[tB]))
        if len(donors) < min_pairs:                      # unlikely, but be safe
            result[f"{tA}‒{tB}"] = np.nan
            continue

        sA = [by_tissue[tA][d] for d in donors]
        sB = [by_tissue[tB][d] for d in donors]

        # Sub-matrices: rows = regions, cols = donors (aligned order!)
        A  = meth[sA].to_numpy(dtype=np.float32)
        B  = meth[sB].to_numpy(dtype=np.float32)

        # Mask NaNs donor-wise per region
        mask = (~np.isnan(A)) & (~np.isnan(B))
        n_ok = mask.sum(axis=1)

        # Pre-allocate with NaN
        col = np.full(len(meth), np.nan, dtype=np.float32)

        # Compute r only where we have enough donors
        valid_rows = np.where(n_ok >= min_pairs)[0]
        for i in valid_rows:
            ai = A[i, mask[i]]
            bi = B[i, mask[i]]
            if ai.size:        # guard against all-NaN after masking
                col[i] = np.corrcoef(ai, bi)[0, 1]

        result[f"{tA}‒{tB}"] = col

    return result

# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meth", required=True, type=Path,
                    help="CSV: rows=regions, cols=SeqCore_ID (samples)")
    ap.add_argument("--map",  required=True, type=Path,
                    help="Tab-delimited map linking SeqCore_ID→GTEx_ID/Tissue")
    ap.add_argument("--min_pairs", default=3, type=int,
                    help="Minimum donors required for a correlation (default=3)")
    ap.add_argument("--out", required=True, type=Path,
                    help="Output CSV: region × tissue-pair matrix")
    args = ap.parse_args()

    meth, smap      = load_inputs(args.meth, args.map)
    by_tissue       = pick_one_sample_per_donor(smap)

    # (If you prefer **averaging duplicates** instead of selecting the first
    #  sample per donor/tissue, replace pick_one_sample_per_donor() with a
    #  routine that takes the column-wise mean of duplicates in `meth`
    #  before correlation.)

    res             = regionwise_corr(meth, by_tissue, args.min_pairs)
    res.to_csv(args.out)
    print(f"✓  wrote {args.out}  ({res.shape[0]:,} regions × {res.shape[1]} tissue-pairs)")

if __name__ == "__main__":
    main()

