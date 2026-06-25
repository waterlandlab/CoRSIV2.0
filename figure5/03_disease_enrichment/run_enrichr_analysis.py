#!/usr/bin/env python3
"""
Enrichr Enrichment Analysis for Top300 TF Gene Lists
======================================================
Runs Enrichr enrichment analysis for each .txt file in the
top300_tf_lists folder and saves results as TSV files.

Output structure:
    enrichr_results/
        <library_name>/
            <list_name>_<library_name>.tsv

Author: auto-generated
Date: 2026-03-03
"""

import json
import os
import time
import requests

# ── Configuration ──────────────────────────────────────────────────────────────

INPUT_DIR = os.path.join(os.path.dirname(__file__), "top300_tf_lists")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "enrichr_results")

# Gene-set libraries to query.  Add or remove as needed.
GENE_SET_LIBRARIES = [
    # ── Pathways & GO ─────────────────────────────────────────
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023",
    "Reactome_2022",
    "MSigDB_Hallmark_2020",
    # ── TF ChIP-Seq ───────────────────────────────────────────
    "ChEA_2022",
    "ENCODE_TF_ChIP-seq_2015",
    "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
    "Transcription_Factor_PPIs",
    # ── Disease databases (Jensen + DisGeNET + others) ────────
    "Jensen_DISEASES",                   # Jensen Lab legacy (~1,811 disease terms)
    "Jensen_DISEASES_Curated_2025",      # Jensen Lab curated, most reliable
    "Jensen_DISEASES_Experimental_2025", # Jensen Lab experimental evidence
    "DisGeNET",                          # DisGeNET (~9,828 disease-gene associations)
    "OMIM_Disease",                      # OMIM disease genes
    "GWAS_Catalog_2023",                 # GWAS Catalog trait associations
    "ClinVar_2019",                      # ClinVar clinical variants
    "Human_Phenotype_Ontology",          # HPO phenotype terms
]

ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"

# Seconds to wait between API calls to be polite to the server
API_DELAY = 0.5

# ── Helpers ─────────────────────────────────────────────────────────────────────

def submit_gene_list(genes: list[str], description: str = "") -> int:
    """Upload a gene list to Enrichr and return the userListId."""
    genes_str = "\n".join(genes)
    payload = {
        "list": (None, genes_str),
        "description": (None, description),
    }
    response = requests.post(f"{ENRICHR_BASE}/addList", files=payload, timeout=60)
    response.raise_for_status()
    data = response.json()
    return data["userListId"]


def fetch_enrichment(user_list_id: int, library: str) -> list:
    """Retrieve enrichment results for a given library."""
    url = f"{ENRICHR_BASE}/enrich"
    params = {"userListId": user_list_id, "backgroundType": library}
    response = requests.get(url, params=params, timeout=60)
    response.raise_for_status()
    data = response.json()
    return data.get(library, [])


def results_to_tsv(results: list, library: str) -> str:
    """
    Convert raw Enrichr results to a TSV string.

    Column order (matches Enrichr export format):
        Rank, Term, P-value, Odds Ratio, Combined Score,
        Overlapping Genes, Adjusted P-value,
        Old P-value, Old Adjusted P-value
    """
    header = "\t".join([
        "Rank", "Term", "P-value", "Odds_Ratio", "Combined_Score",
        "Overlapping_Genes", "Adjusted_P-value",
        "Old_P-value", "Old_Adjusted_P-value",
    ])
    rows = [header]
    for entry in results:
        rank, term, pval, odds, score, genes, adj_pval, old_pval, old_adj_pval = entry
        genes_str = ";".join(genes)
        rows.append(
            f"{rank}\t{term}\t{pval}\t{odds}\t{score}\t"
            f"{genes_str}\t{adj_pval}\t{old_pval}\t{old_adj_pval}"
        )
    return "\n".join(rows)


def read_gene_list(filepath: str) -> list[str]:
    """Read gene symbols from a plain-text file (one per line)."""
    with open(filepath, "r") as fh:
        genes = [line.strip() for line in fh if line.strip()]
    return genes


# ── Main ────────────────────────────────────────────────────────────────────────

def main():
    # Find all .txt files in the input directory
    txt_files = sorted([
        f for f in os.listdir(INPUT_DIR)
        if f.endswith(".txt")
    ])

    if not txt_files:
        print(f"No .txt files found in {INPUT_DIR}")
        return

    print(f"Found {len(txt_files)} gene list file(s) in {INPUT_DIR}")
    print(f"Will query {len(GENE_SET_LIBRARIES)} gene-set libraries\n")

    # Create output directories up-front
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for lib in GENE_SET_LIBRARIES:
        os.makedirs(os.path.join(OUTPUT_DIR, lib), exist_ok=True)

    total = len(txt_files) * len(GENE_SET_LIBRARIES)
    done = 0

    for txt_file in txt_files:
        list_name = os.path.splitext(txt_file)[0]     # e.g. top300_TFs_black_HepG2
        filepath   = os.path.join(INPUT_DIR, txt_file)
        genes      = read_gene_list(filepath)

        print(f"[{list_name}]  {len(genes)} genes")

        # ── Submit gene list ──
        try:
            user_list_id = submit_gene_list(genes, description=list_name)
            print(f"  Submitted → userListId={user_list_id}")
        except Exception as exc:
            print(f"  ERROR submitting gene list: {exc}")
            done += len(GENE_SET_LIBRARIES)
            continue

        time.sleep(API_DELAY)

        # ── Query each library ──
        for library in GENE_SET_LIBRARIES:
            done += 1
            progress = f"[{done}/{total}]"
            out_path = os.path.join(OUTPUT_DIR, library, f"{list_name}__{library}.tsv")

            # Skip if already downloaded (allows re-running without repeating work)
            if os.path.exists(out_path):
                print(f"  {progress} {library} — already exists, skipping.")
                continue

            try:
                results = fetch_enrichment(user_list_id, library)
                tsv_content = results_to_tsv(results, library)

                with open(out_path, "w") as fh:
                    fh.write(tsv_content)

                print(f"  {progress} {library} — {len(results)} terms → {out_path}")

            except Exception as exc:
                print(f"  {progress} {library} — ERROR: {exc}")

            time.sleep(API_DELAY)

        print()  # blank line between files

    print("Done! Results saved to:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
