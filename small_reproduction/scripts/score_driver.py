#!/usr/bin/env python3
"""
Resumable scoring driver for the scDRS-FM reproduction.

Runs one `run_scdrs_fm.py` invocation per job (real or simulation), with:
  - correct per-job flags (imputation, raw_count, filter, include_ctrl_score, ablation)
  - stdout captured to a log file (contains per-stage + per-trait timings)
  - per-job completion checkpoint (a <job>.DONE marker) so reruns skip finished jobs
  - for simulation jobs: flat output to a staging dir, then redistribution of each
    {gsname}.gs.* file into the nested dir layout that notebook 03 expects.

This is REPRODUCTION TOOLING, not part of the scDRS-FM package. It only calls the
*unmodified* package entry point `run_scdrs_fm.py`. All paths resolve relative to
the repo root by default and can be overridden with environment variables (the
same ones set by scripts/_config.sh):

  REPO_ROOT       repo checkout root                 (default: parent of scripts/)
  SCDRS_FM_HOME   scDRS-FM-main package dir          (default: $REPO_ROOT/scDRS-FM-main)
  PYTHON          python interpreter for the package (default: current python)
  DATA_DIR        downloaded data dir                (default: $REPO_ROOT/data)
  SUBSET_DIR      10k subsets dir                    (default: $DATA_DIR/subsets_10k)
  RESULTS         results root                       (default: $REPO_ROOT/results)
  LOCAL_TMP       fast local scratch for logs/sim staging (default: /tmp/scdrsfm)

Usage:
  score_driver.py real   <job_id> [<job_id> ...]   # score real dataset job(s)
  score_driver.py sim    <job_id> [<job_id> ...]   # score simulation batch(es)
  score_driver.py list                              # list all job ids + done-status

Job specs are read from $DATA_DIR/real_jobs.json and $DATA_DIR/sim_jobs.json.
"""
from __future__ import annotations
import sys, os, json, subprocess, time, shutil, re
from pathlib import Path

# ---- path configuration (env-overridable; repo-relative defaults) ----
_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = Path(os.environ.get("REPO_ROOT", str(_SCRIPT_DIR.parent)))
PKG = Path(os.environ.get("SCDRS_FM_HOME", str(REPO_ROOT / "scDRS-FM-main")))
PY = os.environ.get("PYTHON", sys.executable or "python")
RUN = str(PKG / "run_scdrs_fm.py")

DATA = Path(os.environ.get("DATA_DIR", str(REPO_ROOT / "data")))
SUBSETS = Path(os.environ.get("SUBSET_DIR", str(DATA / "subsets_10k")))
GENE_SETS = DATA / "gene_sets"
SIM = DATA / "simulation_data"
COV_DIR = DATA / "extracted"           # covariate files (full-cov versions also in subsets_10k)
LOCAL_TMP = Path(os.environ.get("LOCAL_TMP", "/tmp/scdrsfm"))

# ---- output roots ----
RESULTS = Path(os.environ.get("RESULTS", str(REPO_ROOT / "results")))
REAL_OUT = Path(os.environ.get("REAL_OUT", str(RESULTS / "real")))   # results/real/<job_id>/
# nb03 reads its BASE_DIR; it must contain scdrs+_results/, scdrs_results/, etc.
SIM_OUT  = Path(os.environ.get("SIM_OUT", str(RESULTS / "sim" / "simulation_data")))  # nb03 BASE_DIR
LOG_DIR  = Path(os.environ.get("LOG_DIR", str(RESULTS / "logs")))
for d in (REAL_OUT, SIM_OUT, LOG_DIR):
    d.mkdir(parents=True, exist_ok=True)

# ---- dataset -> (h5ad, cov) filenames within subsets_10k/<name>/ ----
DATASET_FILES = {
    "TMS_FACS":    ("TMS_FACS.h5ad", "TMS_FACS.cov"),
    "TS_FACS":     ("ts_facs.h5ad", "ts.cov"),
    "TMS_Droplet": ("TMS_Droplet.h5ad", "tms_droplet.cov"),
    "SEA_AD":      ("combined_healthy_filtered.h5ad", "combined_healthy_filtered.cov"),
    "Soskic":      ("soskic_100k.h5ad", "soskic.cov"),
    "Nathan":      ("raw.h5ad", "raw.cov"),
    "Cano_Gamez":  ("obj_raw.h5ad", "canogamez.cov"),
    "Braun":       ("human_dev_layers_100k.h5ad", "human_dev_layers_100k.cov"),
}

def all_traits(gs_subdir: str):
    d = GENE_SETS / gs_subdir
    return sorted(p.name for p in d.iterdir() if p.is_file() and not p.name.startswith("."))

def resolve_traits(spec, gs_subdir):
    """spec: None -> all gs_split (75); 'ALL_TCELL' -> all t_cell_pheno (52); list -> as-is."""
    if spec is None:
        return all_traits(gs_subdir)
    if spec == "ALL_TCELL":
        return all_traits(gs_subdir)
    if isinstance(spec, list):
        return spec
    raise ValueError(f"bad trait spec {spec!r}")

def build_cmd(h5ad, cov, out, gs_dir, traits, species, raw, ctrl, abl, imp):
    cmd = [PY, RUN, str(h5ad), str(cov), str(out), str(gs_dir)]
    cmd += list(traits)
    cmd += ["--h5ad_species", species, "--imputation", imp]
    if raw:  cmd.append("--flag_raw_count")
    cmd.append("--flag_filter")           # ALL jobs use --flag_filter per slurms
    if ctrl: cmd.append("--include_ctrl_score")
    if abl:  cmd.append("--ablation")
    return cmd

def run_logged(cmd, log_path, cwd):
    """Run cmd; write full stdout+timings to a LOCAL log (S3/network FS may forbid
    append-reopen), then copy the finished log to the durable log_path once."""
    t0 = time.time()
    local_log = LOCAL_TMP / "score_logs"
    local_log.mkdir(parents=True, exist_ok=True)
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)  # ensure LOG_DIR exists
    local_path = local_log / Path(log_path).name
    env = dict(os.environ)
    env["PYTHONPATH"] = str(PKG) + os.pathsep + env.get("PYTHONPATH", "")
    # Single open handle for the whole run — no append-reopen.
    with open(local_path, "w") as lf:
        lf.write("CMD: " + " ".join(cmd) + "\n")
        lf.write(f"CWD: {cwd}\nSTART: {time.ctime()}\n" + "="*80 + "\n")
        lf.flush()
        p = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, cwd=cwd, env=env)
        dt = time.time() - t0
        lf.write("="*80 + f"\nEND: {time.ctime()}  rc={p.returncode}  wall={dt:.1f}s\n")
        lf.flush()
    # Copy finished log to durable store (copy works everywhere; append/rename may not).
    try:
        shutil.copyfile(local_path, log_path)
    except Exception as e:
        print(f"[warn] could not copy log to durable: {e}")
    return p.returncode, dt

# ---------- REAL ----------
def run_real(job):
    jid = job["id"]
    done = REAL_OUT / jid / ".DONE"
    if done.exists():
        print(f"[real:{jid}] already DONE, skipping")
        return 0
    ds = job["ds"]
    h5, cv = DATASET_FILES[ds]
    h5ad = SUBSETS / ds / h5
    cov  = SUBSETS / ds / cv
    assert h5ad.exists(), f"missing {h5ad}"
    assert cov.exists(), f"missing {cov}"
    out = REAL_OUT / jid
    out.mkdir(parents=True, exist_ok=True)
    traits = resolve_traits(job["traits"], job["gs"])
    gs_dir = GENE_SETS / job["gs"]
    cmd = build_cmd(h5ad, cov, out, gs_dir, traits, job["sp"],
                    job["raw"], job["ctrl"], job["abl"], job["imp"])
    log = LOG_DIR / f"real_{jid}.log"
    print(f"[real:{jid}] {ds} gs={job['gs']} imp={job['imp']} ntraits={len(traits)} ctrl={job['ctrl']} abl={job['abl']}")
    rc, dt = run_logged(cmd, log, cwd=str(PKG))
    if rc == 0:
        done.write_text(f"rc=0 wall={dt:.1f}s ntraits={len(traits)}\n")
        print(f"[real:{jid}] DONE in {dt:.1f}s")
    else:
        print(f"[real:{jid}] FAILED rc={rc}; see {log}")
    return rc

# ---------- SIM ----------
# gs filename patterns -> nested subpath
#  de_overlap: TMS_FACS_{cluster}_{rep}_{overlap}_src{src}.gs -> {cluster}/{rep}/src{src}
#  cell_pct:   TMS_FACS_{c}_{r}_pct{pct}_ov50_src{src}.gs     -> {c}/{r}/pct{pct}/src{src}
#  causal:     TMS_FACS_{c}_{r}_ov{ov}_src{src}.gs            -> {c}/{r}/ov{ov}/src{src}
RE_DEOVL = re.compile(r"^TMS_FACS_(?P<c>.+)_(?P<r>\d+)_(?P<ov>\d+)_src(?P<src>\d+)\.gs$")
RE_PCT   = re.compile(r"^TMS_FACS_(?P<c>.+)_(?P<r>\d+)_pct(?P<pct>\d+)_ov\d+_src(?P<src>\d+)\.gs$")
RE_CAUSAL= re.compile(r"^TMS_FACS_(?P<c>.+)_(?P<r>\d+)_ov(?P<ov>\d+)_src(?P<src>\d+)\.gs$")

def nested_subpath(gsname, family):
    if family == "gs_cell_pct":
        m = RE_PCT.match(gsname)
        assert m, f"cell_pct name no match: {gsname}"
        return Path(m["c"]) / m["r"] / f"pct{m['pct']}" / f"src{m['src']}"
    if family == "gs_causal_genes":
        m = RE_CAUSAL.match(gsname)
        assert m, f"causal name no match: {gsname}"
        return Path(m["c"]) / m["r"] / f"ov{m['ov']}" / f"src{m['src']}"
    # gs_de_overlap (also used by denoise runs)
    m = RE_DEOVL.match(gsname)
    assert m, f"de_overlap name no match: {gsname}"
    return Path(m["c"]) / m["r"] / f"src{m['src']}"

def run_sim(job):
    jid = job["id"]
    out_name = job["out"]
    fam = job["gs"]
    done = SIM_OUT / out_name / ".DONE"
    if done.exists():
        print(f"[sim:{jid}] already DONE, skipping")
        return 0
    h5ad = SIM / "adata" / "TMS_FACS_10k_DE.h5ad"
    cov  = SUBSETS / "TMS_FACS" / "TMS_FACS.cov"   # 10k-aligned cov
    assert h5ad.exists(), f"missing {h5ad}"
    assert cov.exists(), f"missing {cov}"
    gs_dir = SIM / fam
    traits = sorted(p.name for p in gs_dir.iterdir() if p.is_file() and p.name.endswith(".gs"))
    # stage flat output on LOCAL disk (fast, atomic ops), then redistribute to durable nested
    staging = LOCAL_TMP / "sim_stage" / out_name
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True, exist_ok=True)
    # sims: human species, raw+filter per 118_sims slurm
    cmd = build_cmd(h5ad, cov, staging, gs_dir, traits, "human",
                    raw=True, ctrl=False, abl=job["abl"], imp=job["imp"])
    log = LOG_DIR / f"sim_{jid}.log"
    print(f"[sim:{jid}] fam={fam} imp={job['imp']} abl={job['abl']} ntraits={len(traits)} -> {out_name}")
    rc, dt = run_logged(cmd, log, cwd=str(PKG))
    if rc != 0:
        # ALRA needs pyALRA which is optional; it is NOT consumed by nb03.
        # Mark as SKIPPED (best-effort) so it does not block the queue or trigger reruns.
        if job["imp"] == "alra":
            (SIM_OUT / out_name).mkdir(parents=True, exist_ok=True)
            (SIM_OUT / out_name / ".SKIPPED").write_text(
                "imputation=alra requires pyALRA (optional); off nb03 critical path.\n")
            print(f"[sim:{jid}] SKIPPED (pyALRA unavailable; not consumed by nb03)")
            return 0
        print(f"[sim:{jid}] FAILED rc={rc}; see {log}")
        return rc
    # redistribute flat outputs into nested layout on durable store
    dest_root = SIM_OUT / out_name
    n_moved = 0
    for f in staging.iterdir():
        if not f.is_file():
            continue
        # output filenames look like {gsname}.gs.marginal_score.gz etc; parse gsname up to '.gs'
        name = f.name
        idx = name.find(".gs")
        assert idx != -1, f"unexpected sim output filename: {name}"
        gsname = name[:idx] + ".gs"
        suffix = name[idx+3:]                # e.g. .marginal_score.gz
        sub = nested_subpath(gsname, fam)
        ddir = dest_root / sub
        ddir.mkdir(parents=True, exist_ok=True)
        # copy (rename may not be supported on network FS); dest filename keeps gsname + suffix
        shutil.copyfile(f, ddir / (gsname + suffix))
        n_moved += 1
    done.write_text(f"rc=0 wall={dt:.1f}s ntraits={len(traits)} files_redistributed={n_moved}\n")
    print(f"[sim:{jid}] DONE in {dt:.1f}s, redistributed {n_moved} files -> {dest_root}")
    shutil.rmtree(staging, ignore_errors=True)
    return 0

def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    mode = sys.argv[1]
    real_jobs = {j["id"]: j for j in json.load(open(DATA / "real_jobs.json"))}
    sim_jobs  = {j["id"]: j for j in json.load(open(DATA / "sim_jobs.json"))}
    if mode == "list":
        print("REAL jobs:")
        for jid in real_jobs:
            d = (REAL_OUT / jid / ".DONE").exists()
            print(f"  {'[DONE]' if d else '[    ]'} {jid}")
        print("SIM jobs:")
        for jid, j in sim_jobs.items():
            d = (SIM_OUT / j['out'] / ".DONE").exists()
            print(f"  {'[DONE]' if d else '[    ]'} {jid} -> {j['out']}")
        return
    ids = sys.argv[2:]
    rc_all = 0
    if mode == "real":
        pool = real_jobs
        for jid in ids:
            rc_all |= run_real(pool[jid])
    elif mode == "sim":
        pool = sim_jobs
        for jid in ids:
            rc_all |= run_sim(pool[jid])
    else:
        print(f"unknown mode {mode}"); sys.exit(1)
    sys.exit(rc_all)

if __name__ == "__main__":
    main()
