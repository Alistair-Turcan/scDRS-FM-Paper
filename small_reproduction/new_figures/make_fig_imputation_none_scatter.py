#!/usr/bin/env python3
"""
make_fig_imputation_none_scatter.py  --  Figure B of the reproduction:
scDRS-FM (MAGIC imputation) vs standard scDRS (imputation = none) on TMS FACS
(10,000 cells, 75 GWAS traits).

Panel A : pooled per-cell normalized disease score, scDRS-FM (y) vs standard
          scDRS (x), hexbin over all 75 traits x 10k cells, with pooled Pearson r.
Panel B : per-trait detection sensitivity -- % of cells that are FDR-significant
          (P < 0.1) under each method, one point per trait, with the count of
          traits where scDRS-FM >= standard.

Both panels are recomputed from the scored result files shipped in results/:
    scDRS-FM  : $RESULTS/real/tms_facs/*.marginal_score.gz          (MAGIC)
    standard  : $RESULTS/ct/tms_facs_none_ctrl/*.marginal_score.gz  (none)

Outputs (into this script's directory, new_figures/):
    figB_scdrsfm_vs_scdrs_scatter.png / .svg
    figB_scdrsfm_vs_scdrs_per_trait_summary.csv

Paths are env-overridable (RESULTS, or SCDRSFM_BASE as a repo-root fallback);
defaults are repo-relative so a fresh clone with populated results/ just works.

Usage:
    python new_figures/make_fig_imputation_none_scatter.py
"""
from __future__ import annotations
import os, glob
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm  # noqa: F401 (kept for parity with source notebook)

# ---- path resolution (repo-relative defaults; env-overridable) ----
_HERE = Path(__file__).resolve().parent            # new_figures/
_REPO = Path(os.environ.get("SCDRSFM_BASE", str(_HERE.parent)))
RESULTS = Path(os.environ.get("RESULTS", str(_REPO / "results")))
OUT = _HERE

FM_DIR  = Path(os.environ.get("FM_DIR",  str(RESULTS / "real" / "tms_facs")))          # scDRS-FM (MAGIC)
STD_DIR = Path(os.environ.get("STD_DIR", str(RESULTS / "ct" / "tms_facs_none_ctrl")))  # standard scDRS (none)

FDR_P = 0.1  # scDRS FDR-significance threshold used in the paper


def load_pooled():
    fm_traits  = {os.path.basename(f).replace(".marginal_score.gz", "")
                  for f in glob.glob(str(FM_DIR / "*.marginal_score.gz"))}
    std_traits = {os.path.basename(f).replace(".marginal_score.gz", "")
                  for f in glob.glob(str(STD_DIR / "*.marginal_score.gz"))}
    traits = sorted(fm_traits & std_traits)
    if not traits:
        raise SystemExit(
            f"No overlapping *.marginal_score.gz found.\n  FM_DIR = {FM_DIR}\n  STD_DIR = {STD_DIR}\n"
            "Populate results/ (run scripts/01_run_tms_facs.sh and scripts/13_make_ct_tables.sh) "
            "or set FM_DIR/STD_DIR/RESULTS.")
    print(f"traits common to both: {len(traits)}")

    pooled_fm, pooled_std, rows = [], [], []
    for t in traits:
        fm  = pd.read_csv(FM_DIR / f"{t}.marginal_score.gz",  sep="\t", index_col=0)
        std = pd.read_csv(STD_DIR / f"{t}.marginal_score.gz", sep="\t", index_col=0)
        common = fm.index.intersection(std.index)
        a = fm.loc[common, "norm_score"].to_numpy()
        b = std.loc[common, "norm_score"].to_numpy()
        m = np.isfinite(a) & np.isfinite(b)
        a, b = a[m], b[m]
        pooled_fm.append(a); pooled_std.append(b)
        r = np.corrcoef(a, b)[0, 1] if len(a) > 2 else np.nan
        fm_sig  = (fm.loc[common, "pval"]  < FDR_P).mean()
        std_sig = (std.loc[common, "pval"] < FDR_P).mean()
        rows.append({"trait": t, "pearson_r": r,
                     "frac_sig_scdrsfm": fm_sig, "frac_sig_scdrs": std_sig,
                     "n_cells": len(a)})
    pooled_fm  = np.concatenate(pooled_fm)
    pooled_std = np.concatenate(pooled_std)
    summ = pd.DataFrame(rows).sort_values("pearson_r", ascending=False).reset_index(drop=True)
    print("pooled cells:", len(pooled_fm))
    print(f"median per-trait Pearson r: {summ['pearson_r'].median():.3f}")
    print(f"pooled Pearson r: {np.corrcoef(pooled_fm, pooled_std)[0,1]:.3f}")
    return pooled_fm, pooled_std, summ


def make_figure(pooled_fm, pooled_std, summ):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))

    # ---- Panel A: pooled per-cell norm_score hexbin ----
    ax = axes[0]
    lim_lo = np.percentile(np.concatenate([pooled_std, pooled_fm]), 0.2)
    lim_hi = np.percentile(np.concatenate([pooled_std, pooled_fm]), 99.8)
    hb = ax.hexbin(pooled_std, pooled_fm, gridsize=70, bins="log", cmap="viridis",
                   extent=(lim_lo, lim_hi, lim_lo, lim_hi), mincnt=1)
    ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], color="#d62728", lw=1.6, ls="--", label="y = x (identity)")
    cb = fig.colorbar(hb, ax=ax, pad=0.02); cb.set_label("cells per bin (log$_{10}$)")
    r_pool = np.corrcoef(pooled_fm, pooled_std)[0, 1]
    ax.set_xlabel("standard scDRS normalized score\n(imputation = none)")
    ax.set_ylabel("scDRS-FM normalized score\n(MAGIC imputation)")
    ax.set_title(f"A. Per-cell disease scores\n75 traits x 10,000 cells (n={len(pooled_fm):,});  pooled r = {r_pool:.3f}",
                 fontsize=11)
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    ax.set_xlim(lim_lo, lim_hi); ax.set_ylim(lim_lo, lim_hi)
    ax.set_aspect("equal", adjustable="box")

    # ---- Panel B: per-trait significant-cell fraction ----
    ax = axes[1]
    x = summ["frac_sig_scdrs"].to_numpy() * 100
    y = summ["frac_sig_scdrsfm"].to_numpy() * 100
    ax.scatter(x, y, s=26, c="#0279EE", alpha=0.75, edgecolor="white", linewidth=0.4, zorder=3)
    mx = max(x.max(), y.max()) * 1.08
    ax.plot([0, mx], [0, mx], color="#d62728", lw=1.6, ls="--", label="y = x (identity)", zorder=2)
    for _, rr in pd.concat([summ.nlargest(2, "frac_sig_scdrsfm"),
                            summ.nsmallest(1, "frac_sig_scdrsfm")]).iterrows():
        ax.annotate(rr["trait"].replace("UKB_460K.", "").replace("PASS_", ""),
                    (rr["frac_sig_scdrs"] * 100, rr["frac_sig_scdrsfm"] * 100),
                    fontsize=7, xytext=(4, 3), textcoords="offset points", color="#333333")
    n_up = int((y > x).sum())
    n_tot = len(summ)
    ax.set_xlabel("% FDR-significant cells\nstandard scDRS (P < 0.1)")
    ax.set_ylabel("% FDR-significant cells\nscDRS-FM (P < 0.1)")
    ax.set_title(f"B. Detection sensitivity per trait (n={n_tot})\nscDRS-FM >= standard for {n_up}/{n_tot} traits",
                 fontsize=11)
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    ax.set_xlim(0, mx); ax.set_ylim(0, mx)
    ax.set_aspect("equal", adjustable="box")

    fig.suptitle("scDRS-FM vs standard scDRS - TMS FACS (10k cells)", fontsize=13, y=1.02, fontweight="bold")
    fig.tight_layout()
    for ext in ["png", "svg"]:
        fig.savefig(OUT / f"figB_scdrsfm_vs_scdrs_scatter.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    summ.to_csv(OUT / "figB_scdrsfm_vs_scdrs_per_trait_summary.csv", index=False)
    print(f"saved {OUT/'figB_scdrsfm_vs_scdrs_scatter.png'} (+ .svg) and per-trait summary CSV")
    print(f"scDRS-FM >= standard sig-fraction for {n_up}/{n_tot} traits")


if __name__ == "__main__":
    pooled_fm, pooled_std, summ = load_pooled()
    make_figure(pooled_fm, pooled_std, summ)
