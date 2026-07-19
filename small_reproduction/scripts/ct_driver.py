#!/usr/bin/env python
"""
ct_driver.py -- generate .scdrs_ct.<biocol> cell-type association files for
notebooks 01/06/09/10 (full reproduction).

Two phases:
  Phase A (rescore): re-run scDRS-FM WITH --include_ctrl_score for the datasets/
    variants that lack marginal_score_full.gz. Outputs to results/ct/<id>/.
    Variants: magic = scDRS-FM ; none = standard-scDRS baseline.
  Phase B (downstream): for each (result_dir, biocol) run
    scdrs.method.downstream_group_analysis and write <trait><suffix> where
    suffix is ".scdrs_ct.{biocol}" (nb01/06/10) or
    ".marginal_score.gz.scdrs_ct.{biocol}" (nb09).

The 3 datasets already scored with ctrl (nathan_immune, canogamez_immune,
braun_micro) skip Phase A and read full scores from results/real/<id>/.

This is REPRODUCTION TOOLING, not part of the scDRS-FM package. It reuses the
constants/helpers from score_driver.py (same dir) and calls the *unmodified*
package entry point via build_cmd(). Paths are env-overridable exactly as in
score_driver.py / _config.sh; the job matrix is read from
$MATRIX_JSON (default: scripts/scdrs_ct_matrix.json).

Usage:
  ct_driver.py rescore <id>... | rescore-all
  ct_driver.py downstream <result_dir>... | downstream-all
  ct_driver.py list
"""
import sys, os, json, time, shutil, subprocess
from pathlib import Path

# reuse constants + helpers from score_driver (same scripts/ dir)
_SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPT_DIR))
from score_driver import (REPO_ROOT, RESULTS, PKG, PY, RUN, SUBSETS, GENE_SETS,
                          DATASET_FILES, LOG_DIR, resolve_traits, build_cmd,
                          run_logged, all_traits, LOCAL_TMP)

MATRIX_JSON = Path(os.environ.get("MATRIX_JSON", str(_SCRIPT_DIR / "scdrs_ct_matrix.json")))
MATRIX = json.load(open(MATRIX_JSON))
CT_OUT = Path(os.environ.get("CT_OUT", str(RESULTS / "ct")))   # rescored (ctrl) outputs live here
CT_OUT.mkdir(parents=True, exist_ok=True)

# result_dir name -> where its *.marginal_score_full.gz actually live
def resolve_result_dir(name):
    # the 3 "already_have_full" ids live under results/real/<id>/
    already = {j["id"] for j in MATRIX["already_have_full"]}
    if name in already:
        return RESULTS / "real" / name
    return CT_OUT / name

def trait_list_for(spec):
    if spec in (None, "null"):
        return None
    if isinstance(spec, str) and spec in MATRIX["trait_groups"]:
        return list(MATRIX["trait_groups"][spec])
    if isinstance(spec, list):
        return spec
    raise ValueError(f"bad trait spec {spec!r}")

# ---------------- Phase A: rescore with ctrl ----------------
def run_rescore(job):
    jid = job["id"]
    out = CT_OUT / jid
    done = out / ".DONE"
    if done.exists():
        print(f"[rescore:{jid}] already DONE, skipping")
        return 0
    ds = job["ds"]
    h5, cv = DATASET_FILES[ds]
    h5ad = SUBSETS / ds / h5
    cov  = SUBSETS / ds / cv
    assert h5ad.exists(), f"missing {h5ad}"
    assert cov.exists(),  f"missing {cov}"
    out.mkdir(parents=True, exist_ok=True)
    tspec = trait_list_for(job["traits"])
    traits = resolve_traits(None, job["gs"]) if tspec is None else tspec
    gs_dir = GENE_SETS / job["gs"]
    # ctrl ALWAYS True here; abl False
    cmd = build_cmd(h5ad, cov, out, gs_dir, traits, job["sp"],
                    job["raw"], True, False, job["imp"])
    log = LOG_DIR / f"ct_rescore_{jid}.log"
    print(f"[rescore:{jid}] {ds} gs={job['gs']} imp={job['imp']} ntraits={len(traits)} ctrl=True")
    rc, dt = run_logged(cmd, log, cwd=str(PKG))
    if rc == 0:
        # sanity: expect marginal_score_full.gz files
        nfull = len(list(out.glob("*.marginal_score_full.gz")))
        done.write_text(f"rc=0 wall={dt:.1f}s ntraits={len(traits)} nfull={nfull}\n")
        print(f"[rescore:{jid}] DONE in {dt:.1f}s  ({nfull} full-score files)")
    else:
        print(f"[rescore:{jid}] FAILED rc={rc}; see {log}")
    return rc

# ---------------- Phase B: downstream ----------------
_ADATA_CACHE = {}
def load_adata_for_downstream(ds):
    """Load subset h5ad, normalize + log1p + pca + neighbors (scdrs downstream recipe)."""
    if ds in _ADATA_CACHE:
        return _ADATA_CACHE[ds]
    import scanpy as sc, anndata as ad, numpy as np
    h5, _ = DATASET_FILES[ds]
    p = SUBSETS / ds / h5
    print(f"    [adata] loading {p.name} ...")
    a = ad.read_h5ad(p)
    # scDRS downstream: size-factor-normalized + log1p; then pca(20) + neighbors(15,20).
    # Datasets are raw-count subsets.
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.pca(a, n_comps=20)
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=20)
    _ADATA_CACHE[ds] = a
    print(f"    [adata] {ds}: {a.n_obs} cells, connectivities in obsp: {'connectivities' in a.obsp}")
    return a

def run_downstream(dsjob):
    import pandas as pd, scdrs
    result_dir = resolve_result_dir(dsjob["result_dir"])
    ds = dsjob["ds"]
    biocol = dsjob["biocol"]
    suffix = dsjob["suffix"].format(biocol=biocol)
    marker = result_dir / f".DOWNSTREAM_DONE.{biocol}"
    full_files = sorted(result_dir.glob("*.marginal_score_full.gz"))
    if not full_files:
        print(f"[downstream:{dsjob['result_dir']}] NO full-score files in {result_dir} -- did rescore run? SKIP")
        return 1
    if marker.exists():
        print(f"[downstream:{dsjob['result_dir']}] already DONE ({biocol}), skipping")
        return 0
    print(f"[downstream:{dsjob['result_dir']}] ds={ds} biocol={biocol} nfull={len(full_files)} suffix='{suffix}'")
    adata = load_adata_for_downstream(ds)
    assert biocol in adata.obs.columns, f"{biocol} not in {ds}.obs"
    n_ok = 0
    for f in full_files:
        prefix = f.name.replace(".marginal_score_full.gz", "")
        out_path = result_dir / f"{prefix}{suffix}"
        if out_path.exists():
            n_ok += 1; continue
        df_full = pd.read_csv(f, sep="\t", index_col=0)
        try:
            res = scdrs.method.downstream_group_analysis(
                adata=adata, df_full_score=df_full, group_cols=[biocol])
        except Exception as e:
            print(f"    [FAIL] {prefix}: {e}")
            continue
        dfres = res[biocol]
        # write with cell-type as a column named exactly biocol (nb readers do read_csv sep=\t)
        dfres = dfres.copy()
        dfres.index.name = biocol
        # write to LOCAL then copy (network FS safe for text; .scdrs_ct is plain tsv)
        tmp = LOCAL_TMP / "ct_tmp"; tmp.mkdir(parents=True, exist_ok=True)
        lp = tmp / out_path.name
        dfres.to_csv(lp, sep="\t")
        shutil.copyfile(lp, out_path)
        n_ok += 1
    marker.write_text(f"biocol={biocol} nfiles={n_ok} suffix={suffix}\n")
    print(f"[downstream:{dsjob['result_dir']}] wrote {n_ok}/{len(full_files)} .scdrs_ct files")
    return 0

def cmd_list():
    print("=== RESCORE JOBS ===")
    for j in MATRIX["rescore_jobs"]:
        out = CT_OUT / j["id"]
        st = "DONE" if (out/".DONE").exists() else "----"
        nfull = len(list(out.glob("*.marginal_score_full.gz"))) if out.exists() else 0
        print(f"  [{st}] {j['id']:<26} {j['ds']:<11} imp={j['imp']:<5} full={nfull}")
    print("=== DOWNSTREAM JOBS ===")
    for j in MATRIX["downstream_jobs"]:
        rd = resolve_result_dir(j["result_dir"])
        bio = j["biocol"]
        st = "DONE" if (rd/f".DOWNSTREAM_DONE.{bio}").exists() else "----"
        suffix = j["suffix"].format(biocol=bio)
        nct = len(list(rd.glob(f"*{suffix}"))) if rd.exists() else 0
        print(f"  [{st}] {j['result_dir']:<26} {j['ds']:<11} biocol={bio:<18} nct={nct}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    mode = sys.argv[1]
    rc = 0
    if mode == "rescore":
        jobs = {j["id"]: j for j in MATRIX["rescore_jobs"]}
        for jid in sys.argv[2:]:
            rc |= run_rescore(jobs[jid])
    elif mode == "rescore-all":
        for j in MATRIX["rescore_jobs"]:
            rc |= run_rescore(j)
    elif mode == "downstream":
        jobs = {j["result_dir"]: j for j in MATRIX["downstream_jobs"]}
        for name in sys.argv[2:]:
            rc |= run_downstream(jobs[name])
    elif mode == "downstream-all":
        for j in MATRIX["downstream_jobs"]:
            rc |= run_downstream(j)
    elif mode == "list":
        cmd_list()
    else:
        print(__doc__); sys.exit(1)
    sys.exit(rc)
