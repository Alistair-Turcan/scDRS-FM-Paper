#!/usr/bin/env python3
"""
make_fig_runtime_memory.py  --  Figure A of the reproduction: scDRS-FM
computational cost on TMS FACS (10,000 cells, 75 GWAS traits, MAGIC imputation).

Four panels:
  A. One-time preprocessing phases (horizontal bars): scDRS preprocess, metacell
     assignment, MAGIC imputation, covariate correction, metacell aggregation.
     These run ONCE per dataset, independent of the number of traits.
  B. Per-trait scoring-time histogram (n = 75 traits) with mean + median.
  C. Total scoring budget (one-time vs per-trait scoring, stacked) + peak RSS box.
  D. Mean per-trait scoring-time composition (build gene set + control, marginal,
     conditional, save marginal, variance-ratio).

DATA SOURCES (both ship in new_figures/):
  * per-trait timings : figA_per_trait_timings.csv
        columns: trait, build_gs_ctrl, var_ratio, marginal, save_marginal,
                 conditional, trait_total
        Parsed from the definitive scoring log's per-trait timing table
        (scripts/01_run_tms_facs.sh -> results/logs/real_tms_facs.log).
  * peak memory       : memprofile_result.json  (key: peak_rss_gb)
        A REAL measurement from scripts/memprofile_scoring.py using
        resource.getrusage(RUSAGE_CHILDREN).ru_maxrss on a 2-trait scoring run.

ONE-TIME PHASE TIMINGS are constants below, parsed from the same definitive log
(real_tms_facs.log). To re-derive them from a log, set SCDRSFM_SCORING_LOG to a
log path and re-run; otherwise the documented constants are used.

Outputs (into this script's directory, new_figures/):
    figA_scdrsfm_runtime_memory.png / .svg

Usage:
    python new_figures/make_fig_runtime_memory.py
"""
from __future__ import annotations
import os, re, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent
TIMINGS_CSV = Path(os.environ.get("FIGA_TIMINGS_CSV", str(_HERE / "figA_per_trait_timings.csv")))
MEMPROFILE  = Path(os.environ.get("FIGA_MEMPROFILE", str(_HERE / "memprofile_result.json")))
SCORING_LOG = os.environ.get("SCDRSFM_SCORING_LOG")  # optional: re-derive one-time phases

# ---- One-time preprocessing phase timings (seconds), from the definitive
#      real_tms_facs.log (scripts/01_run_tms_facs.sh). Documented constants. ----
ONETIME_DEFAULT = {
    "scDRS preprocess":     17.846,
    "Metacell assignment":  65.425,
    "MAGIC imputation":     58.061,
    "Covariate correction":  2.403,
    "Metacell aggregation":  1.331,
}

# Colorblind-friendly palette (Phylo)
C_MAGIC  = "#0279EE"
C_META   = "#75A025"
C_PREP   = "#FF9400"
C_SCORE  = "#7A7A7A"
C_ACCENT = "#D45E00"


def parse_onetime_from_log(logpath: str) -> dict:
    """Best-effort re-derivation of one-time phase totals from a scoring log.
    Falls back to ONETIME_DEFAULT for any phase not found."""
    text = Path(logpath).read_text().splitlines()

    def find_val(patterns):
        for pat in patterns:
            for line in text:
                m = re.search(pat, line)
                if m:
                    return float(m.group(1))
        return None

    got = dict(ONETIME_DEFAULT)
    mapping = {
        "scDRS preprocess":     [r"[Pp]reprocess(?:ing)?[^\d]*([\d.]+)\s*s"],
        "Metacell assignment":  [r"[Mm]etacell assign[^\d]*([\d.]+)\s*s"],
        "MAGIC imputation":     [r"Calculated MAGIC in ([\d.]+) seconds",
                                 r"[Ii]mputation[^\d]*([\d.]+)\s*s"],
        "Covariate correction": [r"[Cc]ovariate[^\d]*([\d.]+)\s*s"],
        "Metacell aggregation": [r"[Mm]etacell aggregat[^\d]*([\d.]+)\s*s"],
    }
    for k, pats in mapping.items():
        v = find_val(pats)
        if v is not None:
            got[k] = v
    return got


def load_data():
    if not TIMINGS_CSV.exists():
        raise SystemExit(f"Missing per-trait timings CSV: {TIMINGS_CSV}")
    df = pd.read_csv(TIMINGS_CSV)
    onetime = parse_onetime_from_log(SCORING_LOG) if SCORING_LOG else dict(ONETIME_DEFAULT)
    if MEMPROFILE.exists():
        peak_gb = float(json.loads(MEMPROFILE.read_text())["peak_rss_gb"])
    else:
        peak_gb = float("nan")
        print(f"[warn] {MEMPROFILE} not found; peak-RSS box will be omitted.")
    return df, onetime, peak_gb


def make_figure(df, onetime, peak_gb):
    onetime_total = sum(onetime.values())
    n_traits = len(df)

    fig = plt.figure(figsize=(15.5, 8.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.42], width_ratios=[1.0, 1.05, 0.85],
                          hspace=0.42, wspace=0.40)

    # ===== Panel A: one-time phases (horizontal bars) =====
    axA = fig.add_subplot(gs[0, 0])
    phase_order = ["scDRS preprocess", "Metacell assignment", "MAGIC imputation",
                   "Covariate correction", "Metacell aggregation"]
    phase_vals = [onetime[p] for p in phase_order]
    phase_colors = [C_PREP, C_META, C_MAGIC, "#FFB84D", "#A6C97A"]
    ypos = np.arange(len(phase_order))[::-1]
    axA.barh(ypos, phase_vals, color=phase_colors, edgecolor="white", height=0.68)
    for y, v in zip(ypos, phase_vals):
        axA.text(v + 1.2, y, f"{v:.1f}s", va="center", ha="left", fontsize=9.5, fontweight="bold")
    axA.set_yticks(ypos); axA.set_yticklabels(phase_order, fontsize=9.5)
    axA.set_xlabel("Time (seconds)", fontsize=11)
    axA.set_xlim(0, max(phase_vals) * 1.28)
    axA.set_title(f"A. One-time preprocessing phases\n(total {onetime_total:.0f}s, independent of #traits)",
                  fontsize=11, fontweight="bold")
    axA.spines[["top", "right"]].set_visible(False)

    # ===== Panel B: per-trait scoring histogram =====
    axB = fig.add_subplot(gs[0, 1])
    axB.hist(df["trait_total"], bins=18, color=C_SCORE, edgecolor="white", alpha=0.88)
    axB.axvline(df["trait_total"].mean(), color=C_ACCENT, lw=2.2, ls="--",
                label=f"mean = {df['trait_total'].mean():.1f}s")
    axB.axvline(df["trait_total"].median(), color="#0279EE", lw=1.8, ls=":",
                label=f"median = {df['trait_total'].median():.1f}s")
    axB.set_xlabel("Per-trait scoring time (seconds)", fontsize=11)
    axB.set_ylabel("Number of GWAS traits", fontsize=11)
    axB.set_title(f"B. Per-trait scoring time\n(n={n_traits} traits, 10k cells)",
                  fontsize=11, fontweight="bold")
    axB.legend(fontsize=9.5, frameon=False, loc="upper right")
    axB.spines[["top", "right"]].set_visible(False)
    axB.set_xlim(df["trait_total"].min() - 1, df["trait_total"].max() + 1.2)

    # ===== Panel C: total budget + memory =====
    axC = fig.add_subplot(gs[0, 2])
    total_scoring = df["trait_total"].sum()
    axC.bar(0, onetime_total, color=C_MAGIC, edgecolor="white", width=0.55, label="One-time preprocessing")
    axC.bar(0, total_scoring, bottom=onetime_total, color=C_SCORE, edgecolor="white", width=0.55,
            label=f"Scoring ({n_traits} traits)")
    grand = onetime_total + total_scoring
    axC.text(0, onetime_total / 2, f"{onetime_total:.0f}s\n({100*onetime_total/grand:.0f}%)",
             ha="center", va="center", color="white", fontsize=9, fontweight="bold")
    axC.text(0, onetime_total + total_scoring / 2, f"{total_scoring:.0f}s\n({100*total_scoring/grand:.0f}%)",
             ha="center", va="center", color="white", fontsize=9, fontweight="bold")
    axC.text(0, grand + 18, f"total {grand/60:.1f} min", ha="center", va="bottom", fontsize=10, fontweight="bold")
    axC.set_xlim(-0.55, 1.0); axC.set_xticks([])
    axC.set_ylabel("Wall-clock time (seconds)", fontsize=11)
    axC.set_ylim(0, grand * 1.20)
    axC.set_title("C. Total scoring budget\n& peak memory", fontsize=11, fontweight="bold")
    axC.legend(fontsize=8.5, frameon=False, loc="upper right", bbox_to_anchor=(1.05, 0.66))
    axC.spines[["top", "right"]].set_visible(False)
    if np.isfinite(peak_gb):
        axC.text(0.66, grand * 0.30, f"Peak RSS\n{peak_gb:.1f} GB",
                 ha="center", va="center", fontsize=11, fontweight="bold", color=C_ACCENT,
                 bbox=dict(boxstyle="round,pad=0.4", fc="#FFF3E0", ec=C_ACCENT, lw=1.5))

    # ===== Panel D (spans bottom): mean per-trait component composition =====
    axD = fig.add_subplot(gs[1, :])
    comp_cols = ["build_gs_ctrl", "marginal", "conditional", "save_marginal", "var_ratio"]
    comp_labels = ["Build gene set + control genes", "Marginal scoring", "Conditional scoring",
                   "Save marginal", "Variance-ratio"]
    comp_colors = ["#0279EE", "#75A025", "#FF9400", "#B79CED", "#9A9A9A"]
    comp_means = [df[c].mean() for c in comp_cols]
    left = 0
    for lab, v, c in zip(comp_labels, comp_means, comp_colors):
        axD.barh(0, v, left=left, color=c, edgecolor="white", height=0.5, label=f"{lab}: {v:.2f}s")
        if v > 0.6:
            axD.text(left + v / 2, 0, f"{v:.1f}s", ha="center", va="center", color="white",
                     fontsize=9, fontweight="bold")
        left += v
    axD.set_xlim(0, left * 1.02); axD.set_ylim(-0.5, 0.7); axD.set_yticks([])
    axD.set_xlabel("Mean time per trait (seconds)", fontsize=10.5)
    axD.set_title(f"D. Average per-trait scoring time composition (mean total = {df['trait_total'].mean():.1f}s / trait)",
                  fontsize=11, fontweight="bold")
    axD.spines[["top", "right", "left"]].set_visible(False)
    axD.legend(fontsize=9, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.35), ncol=5,
               handlelength=1.2, columnspacing=1.5)

    fig.suptitle("scDRS-FM computational cost - TMS FACS, 10,000 cells, 75 GWAS traits (MAGIC imputation)",
                 fontsize=13, fontweight="bold", y=0.98)

    for ext in ["png", "svg"]:
        fig.savefig(_HERE / f"figA_scdrsfm_runtime_memory.{ext}",
                    dpi=150 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {_HERE/'figA_scdrsfm_runtime_memory.png'} (+ .svg)")
    print(f"one-time={onetime_total:.1f}s  per-trait mean={df['trait_total'].mean():.2f}s x {n_traits}"
          f"  grand total={grand:.1f}s ({grand/60:.1f} min)  peak RSS={peak_gb:.3f} GB")


if __name__ == "__main__":
    df, onetime, peak_gb = load_data()
    make_figure(df, onetime, peak_gb)
