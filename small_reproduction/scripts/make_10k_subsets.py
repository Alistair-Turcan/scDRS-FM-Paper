#!/usr/bin/env python3
"""
Subsample each real scDRS-FM dataset to EXACTLY 10,000 cells (fixed seed=0).

Design rationale (faithful to the scDRS-FM pipeline):
  The pipeline (run_scdrs_fm.py, --flag_filter) applies sc.pp.filter_cells(min_genes=250)
  BEFORE scoring. To guarantee exactly 10k cells reach MAGIC + scoring (not fewer after
  in-pipeline filtering), we subsample from the pool of cells that already pass min_genes>=250.
  We PRESERVE the raw X matrix and the full gene set (gene filtering happens in-pipeline on
  the 10k subset, consistent with the original per-dataset run). The covariate file is
  subset to the same barcodes and written alongside.

Outputs (per dataset) -> $SUBSET_DIR/<DatasetName>/  (default: <repo>/data/subsets_10k/<DatasetName>/)
  - <h5ad basename>            (10k-cell h5ad, raw X preserved)
  - <cov basename>             (10k-row cov, aligned to subset barcodes)
Braun reconciliation: figshare ships human_dev_layers_100k.h5ad + .cov, but the slurm
  references *_gene_symbols.h5ad + *_processed.cov. We write BOTH names (copy) so either
  path resolves.

Uses backed='r' read to avoid loading the full (up to 19.5 GB) matrix into RAM.
"""
import gc
import os
import shutil
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

SEED = 0
TARGET_N = 10_000
MIN_GENES = 250  # matches pipeline filter_cells(min_genes=250)

EXT = os.environ.get("DATA_EXTRACTED", os.path.join(os.environ.get("DATA_DIR", "./data"), "extracted"))
OUT_ROOT = os.environ.get("SUBSET_DIR", os.path.join(os.environ.get("DATA_DIR", "./data"), "subsets_10k"))
LOCAL_TMP = os.environ.get("LOCAL_TMP", "/tmp/subset_tmp")  # h5ad writes need local disk (random-access)

# name -> (h5ad, cov, species, is_raw_count)
DATASETS = [
    ("TMS_FACS",   "TMS_FACS.h5ad",                  "TMS_FACS.cov",                  "mouse", True),
    ("TS_FACS",    "ts_facs.h5ad",                   "ts.cov",                        "human", True),
    ("TMS_Droplet","TMS_Droplet.h5ad",               "tms_droplet.cov",               "mouse", True),
    ("SEA_AD",     "combined_healthy_filtered.h5ad", "combined_healthy_filtered.cov", "human", False),
    ("Soskic",     "soskic_100k.h5ad",               "soskic.cov",                    "human", True),
    ("Nathan",     "raw.h5ad",                       "raw.cov",                       "human", True),
    ("Cano_Gamez", "obj_raw.h5ad",                   "canogamez.cov",                 "human", True),
    ("Braun",      "human_dev_layers_100k.h5ad",     "human_dev_layers_100k.cov",     "human", True),
]


def compute_n_genes_backed(A, chunk=20000):
    """Number of detected genes per cell (nnz per row), chunked to bound memory."""
    n = A.shape[0]
    ng = np.zeros(n, dtype=np.int32)
    for start in range(0, n, chunk):
        end = min(start + chunk, n)
        Xc = A.X[start:end]
        if sp.issparse(Xc):
            ng[start:end] = np.asarray((Xc > 0).sum(axis=1)).ravel()
        else:
            ng[start:end] = (np.asarray(Xc) > 0).sum(axis=1)
    return ng


def main():
    only = sys.argv[1] if len(sys.argv) > 1 else None
    os.makedirs(OUT_ROOT, exist_ok=True)
    os.makedirs(LOCAL_TMP, exist_ok=True)
    summary = []

    for name, h5, cov, species, is_raw in DATASETS:
        if only and name != only:
            continue
        print(f"\n{'='*70}\n{name}: {h5}  (species={species}, raw_count={is_raw})", flush=True)
        h5_path = os.path.join(EXT, h5)
        cov_path = os.path.join(EXT, cov)
        out_dir = os.path.join(OUT_ROOT, name)
        os.makedirs(out_dir, exist_ok=True)

        # --- read backed, compute eligible pool ---
        A = ad.read_h5ad(h5_path, backed="r")
        n_total = A.shape[0]
        print(f"  total cells: {n_total:,}  genes: {A.shape[1]:,}", flush=True)

        # use cov n_genes if present & aligned, else compute from X
        ng = None
        try:
            dfc_head = pd.read_csv(cov_path, sep="\t", index_col=0, nrows=5)
            if "n_genes" in dfc_head.columns:
                dfc_full = pd.read_csv(cov_path, sep="\t", index_col=0)
                # align to adata obs_names
                common = A.obs_names.isin(dfc_full.index)
                if common.all() and len(dfc_full) >= n_total:
                    ng = dfc_full.reindex(A.obs_names)["n_genes"].values.astype(float)
                    if np.isnan(ng).any():
                        ng = None
        except Exception as e:
            print(f"  (cov n_genes read failed: {e})", flush=True)
        if ng is None:
            print("  computing n_genes from X (chunked)...", flush=True)
            ng = compute_n_genes_backed(A)

        eligible = np.where(ng >= MIN_GENES)[0]
        print(f"  eligible (n_genes>={MIN_GENES}): {len(eligible):,}", flush=True)
        if len(eligible) < TARGET_N:
            print(f"  WARNING: only {len(eligible)} eligible < {TARGET_N}; using all eligible", flush=True)
            sel = eligible
        else:
            rng = np.random.default_rng(SEED)
            sel = np.sort(rng.choice(eligible, size=TARGET_N, replace=False))
        print(f"  selected: {len(sel):,} cells", flush=True)

        # --- materialize subset (load only selected rows) ---
        sub = A[sel].to_memory()
        A.file.close()
        del A
        gc.collect()
        # ensure csr + keep raw X as-is
        if sp.issparse(sub.X):
            sub.X = sub.X.tocsr()
        print(f"  subset shape: {sub.shape}", flush=True)

        # write h5ad to LOCAL disk first (random-access), then copy to shared
        local_h5 = os.path.join(LOCAL_TMP, h5)
        sub.write_h5ad(local_h5)
        shutil.copy(local_h5, os.path.join(out_dir, h5))
        os.remove(local_h5)
        print(f"  wrote h5ad -> {out_dir}/{h5}", flush=True)

        # --- subset cov (aligned to subset barcodes) ---
        dfc = pd.read_csv(cov_path, sep="\t", index_col=0)
        dfc_sub = dfc.reindex(sub.obs_names)
        n_missing = dfc_sub.isnull().any(axis=1).sum()
        if n_missing:
            print(f"  WARNING: {n_missing} subset cells missing in cov (barcode mismatch)", flush=True)
        dfc_sub.to_csv(os.path.join(out_dir, cov), sep="\t")
        print(f"  wrote cov  -> {out_dir}/{cov}  ({len(dfc_sub)} rows)", flush=True)

        # --- Braun name reconciliation: also write slurm-expected names ---
        if name == "Braun":
            shutil.copy(os.path.join(out_dir, h5),
                        os.path.join(out_dir, "human_dev_layers_100k_gene_symbols.h5ad"))
            shutil.copy(os.path.join(out_dir, cov),
                        os.path.join(out_dir, "human_dev_layers_100k_processed.cov"))
            print("  Braun: also wrote *_gene_symbols.h5ad + *_processed.cov", flush=True)

        summary.append({
            "name": name, "h5ad": h5, "cov": cov, "species": species,
            "raw_count": is_raw, "n_total": int(n_total),
            "n_eligible": int(len(eligible)), "n_selected": int(len(sel)),
            "n_genes_full": int(sub.shape[1]),
        })
        del sub
        gc.collect()

    # write summary
    import json
    sm_path = os.path.join(OUT_ROOT, "subset_summary.json")
    json.dump({"seed": SEED, "target_n": TARGET_N, "min_genes": MIN_GENES,
               "datasets": summary}, open(sm_path, "w"), indent=2)
    print(f"\n{'='*70}\nSUMMARY written -> {sm_path}", flush=True)
    for s in summary:
        print(f"  {s['name']:12s} {s['n_total']:>8,} -> {s['n_selected']:>6,} cells "
              f"(eligible {s['n_eligible']:,})", flush=True)


if __name__ == "__main__":
    main()
