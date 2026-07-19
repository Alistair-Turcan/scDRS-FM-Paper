#!/usr/bin/env python
"""Memory-profiled wrapper around run_scdrs_fm.py to capture REAL peak RSS.

Runs the identical scDRS-FM scoring pipeline as a subprocess and records the
child's peak resident-set-size via resource.getrusage(RUSAGE_CHILDREN). The
one-time MAGIC/preprocessing phase dominates memory and is independent of trait
count, so a 2-trait run yields a faithful peak-memory estimate for the full
75-trait job. This produces the peak_rss_gb used in Figure A (Panel C).

This is REPRODUCTION TOOLING: it only *calls* the unmodified scDRS-FM package
entry point run_scdrs_fm.py; it does not modify any package code.

Paths are env-overridable exactly as in scripts/_config.sh (repo-relative
defaults): REPO_ROOT, SCDRS_FM_HOME, PYTHON, DATA_DIR, SUBSET_DIR, RESULTS,
LOCAL_TMP. Result JSON is written to $RESULTS/logs/memprofile_result.json and
should be copied next to the figure scripts (new_figures/) for Figure A.

Usage:
    python scripts/memprofile_scoring.py
"""
import os
import json
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = Path(os.environ.get("REPO_ROOT", str(_SCRIPT_DIR.parent)))
PKG = Path(os.environ.get("SCDRS_FM_HOME", str(REPO_ROOT / "scDRS-FM-main")))
PY = os.environ.get("PYTHON", sys.executable or "python")
DATA = Path(os.environ.get("DATA_DIR", str(REPO_ROOT / "data")))
SUBSETS = Path(os.environ.get("SUBSET_DIR", str(DATA / "subsets_10k")))
RESULTS = Path(os.environ.get("RESULTS", str(REPO_ROOT / "results")))
LOCAL_TMP = Path(os.environ.get("LOCAL_TMP", "/tmp/scdrsfm"))

SCRIPT = str(PKG / "run_scdrs_fm.py")
H5AD = str(SUBSETS / "TMS_FACS" / "TMS_FACS.h5ad")
COV = str(SUBSETS / "TMS_FACS" / "TMS_FACS.cov")
GS = str(DATA / "gene_sets" / "gs_split")
OUT = str(LOCAL_TMP / "memprofile_out")   # local scratch, discarded
TRAITS = ["PASS_ADHD_Demontis2018", "PASS_MDD_Howard2019"]

Path(OUT).mkdir(parents=True, exist_ok=True)
(RESULTS / "logs").mkdir(parents=True, exist_ok=True)

cmd = [PY, SCRIPT, H5AD, COV, OUT, GS] + TRAITS + [
    "--h5ad_species", "mouse",
    "--imputation", "magic",
    "--flag_raw_count",
    "--flag_filter",
]

env = dict(os.environ)
env["PYTHONPATH"] = str(PKG) + os.pathsep + env.get("PYTHONPATH", "")

print("CMD:", " ".join(cmd), flush=True)
t0 = time.time()
rc = subprocess.call(cmd, env=env)
wall = time.time() - t0

usage = resource.getrusage(resource.RUSAGE_CHILDREN)
# ru_maxrss units: kilobytes on Linux, bytes on macOS
maxrss_kb = usage.ru_maxrss
if platform.system() == "Darwin":
    maxrss_kb = maxrss_kb / 1024.0
peak_gb = maxrss_kb / (1024.0 ** 2)

result = {
    "rc": rc,
    "wall_s": round(wall, 2),
    "n_traits": len(TRAITS),
    "peak_rss_gb": round(peak_gb, 3),
    "peak_rss_kb": maxrss_kb,
    "n_cells": 10000,
    "note": "one-time MAGIC/preprocessing dominates peak memory; trait-count independent",
}
print("MEMPROFILE_RESULT:", json.dumps(result), flush=True)
out_json = RESULTS / "logs" / "memprofile_result.json"
with open(out_json, "w") as f:
    json.dump(result, f, indent=2)
print(f"Saved to {out_json}", flush=True)
print("Copy this file to new_figures/memprofile_result.json for Figure A.", flush=True)
