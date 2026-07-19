# scDRS-FM Reproduction (10k-cell scale)

A GitHub-ready reproduction of the **scDRS-FM** paper. This repository does **not**
modify the scDRS-FM method: it vendors the released package **unchanged** under
[`scDRS-FM-main/`](scDRS-FM-main/) and drives it with a set of thin, fully
parameterized wrapper scripts and the paper's analysis notebooks. Every result,
figure, and cell-type table in this repo is produced by calling the *unmodified*
`scDRS-FM-main/run_scdrs_fm.py` (and the package's other entry points) exactly as
the original SLURM jobs did.

> **Scale note.** To make the full pipeline runnable end-to-end on a single
> workstation, **every real dataset is subsampled to exactly 10,000 cells**
> (fixed `seed=0`, keeping cells with `>=250` genes and preserving raw counts and
> the full gene set). This is the only substantive deviation from the paper, and
> it is confined to a repo-tooling script (`scripts/make_10k_subsets.py`) — the
> scDRS-FM package itself is untouched. See [`docs/DATA.md`](docs/DATA.md).

---

## 1. What's in here

```
scDRS-FM_reproduction/
├── scDRS-FM-main/        # the scDRS-FM package, VENDORED UNMODIFIED (do not edit)
├── download/             # data + MAGMA reference download scripts + manifest
├── scripts/              # parameterized bash/python drivers (00–13) that CALL the package
├── notebooks/            # the paper's analysis notebooks (executed, WITH outputs)
├── new_figures/          # two new figures added by this reproduction (+ their generators)
├── results/              # figures/ + tables/ shipped; score dirs regenerated (see below)
├── docs/                 # DATA.md (subsampling) + NB_EDIT_MAP.md (per-notebook path edits)
├── requirements.txt      # pip dependencies for the analysis environment
├── environment.yml       # optional conda environment
└── LICENSE
```

* **`scDRS-FM-main/`** — the released scDRS-FM package, byte-for-byte unmodified.
  All scoring goes through its `run_scdrs_fm.py` CLI; all gradient decomposition
  through its `run_decompose_gradients.py`. **Do not edit anything under this
  directory** — that is the whole point of the reproduction.
* **`scripts/`** — orchestration only. These reproduce the original SLURM array
  jobs as plain, resumable bash/python. They construct the *exact* CLI the paper
  used and never reimplement any method logic. All paths default to repo-relative
  locations and are overridable by environment variables (see `scripts/_config.sh`).
* **`notebooks/`** — the paper's figure/analysis notebooks, **shipped executed with
  their outputs** so you can read the reproduced results without rerunning. The
  only edits are input-path retargeting to this repo's layout (documented in
  [`docs/NB_EDIT_MAP.md`](docs/NB_EDIT_MAP.md)); no analysis logic was changed.
* **`new_figures/`** — two figures this reproduction adds (runtime/memory profile
  and a scDRS-FM-vs-scDRS concordance scatter) plus the standalone scripts that
  regenerate them from the score files.

**What `results/` ships vs. regenerates.** The raw score tree (all
`*.marginal_score.gz` etc.) is ~32 GB and is **regenerated** by the scripts, so
`results/real/`, `results/ct/`, `results/sim/`, `results/null/`, and
`results/logs/` ship empty (a `.gitkeep` documents each). What *is* shipped:
`results/figures/` (the two new figures) and `results/tables/` (the manuscript
supplementary cell-type-proportion CSVs). Every notebook figure/table is also
embedded inline in the executed notebooks under `notebooks/`.

---

## 2. Environment setup

Python 3.9–3.11 is recommended. Using a fresh virtual environment:

```bash
# pip
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# make the unmodified package importable
pip install -e scDRS-FM-main            # or: export PYTHONPATH=$PWD/scDRS-FM-main

# --- or conda ---
conda env create -f environment.yml
conda activate scdrsfm-repro
pip install -e scDRS-FM-main
```

Key dependencies (pinned where a version matters): `scdrs`, `magic-impute`,
`scprep`, `scanpy`, `anndata`, `leidenalg`, `gseapy==1.3.0`,
`gprofiler-official`, plus the usual `numpy/scipy/pandas/matplotlib/statsmodels`.
`gseapy==1.3.0` is required by nb09 and nb10; `magic-impute` (with
`random_state=0`) provides the MAGIC imputation that defines scDRS-FM.

**MAGMA** is a separate C++ binary (not a Python package). It is downloaded by
`download/download_magma_reference.sh` and is only needed for step P6 (regenerating
the gene-level Z matrix for nb08).

---

## 3. Get the data

Two downloads. Both write into the repo by default.

```bash
# (a) scDRS-FM data release (figshare article 33000602):
#     8 datasets (h5ad + cov), 75 GWAS gene sets, microglia + T-cell signatures,
#     genet_cor.csv, and 75 LDSC sumstats. ~10 GB compressed.
bash download/download_figshare_data.sh            # -> ./data/{raw,extracted,sumstats,gene_sets}

# (b) MAGMA v1.10 binary + 1000G Phase 3 EUR reference + NCBI37.3 gene locations.
bash download/download_magma_reference.sh          # -> ./magma_ref/
```

The expected file inventory is in [`download/manifest.tsv`](download/manifest.tsv).
If figshare's layout changes, download the files manually and place them under
`data/` following the manifest.

---

## 4. Reproduce, step by step

All scripts source `scripts/_config.sh`, which resolves every path relative to the
repo root (overridable via env vars) and defines the `run_scdrs_fm` bash helper
that shells out to the unmodified package. Run them **from the repo root**.

### P0 — Subsample every dataset to 10k cells
```bash
bash scripts/00_make_10k_subsets.sh                # all datasets -> data/subsets_10k/
```

### P1 — Score the real datasets (scDRS-FM)
Each script reproduces the exact per-dataset flags from the paper's SLURM jobs
(species, filtering, raw-count handling, imputation, and whether control scores
are included). Run the ones you need:

```bash
bash scripts/01_run_tms_facs.sh      # TMS_FACS  (mouse, MAGIC, 75 GWAS traits)  ← headline run
bash scripts/02_run_ts_facs.sh       # TS_FACS   (human, imputation none, +ctrl)
bash scripts/03_run_tms_droplet.sh   # TMS_Droplet (mouse, none, +ctrl)
bash scripts/04_run_soskic.sh        # Soskic    (human, MAGIC): 15 immune GWAS + 52 T-cell sigs
bash scripts/05_run_sea_ad.sh        # SEA_AD    (human, MAGIC, 8 brain GWAS + microglia)
bash scripts/06_run_braun.sh         # Braun     (human, MAGIC +ctrl): 8 brain GWAS + 5 microglia
bash scripts/07_run_nathan.sh        # Nathan    (human, MAGIC +ctrl): immune GWAS + T-cell sigs
bash scripts/08_run_canogamez.sh     # Cano_Gamez(human, MAGIC +ctrl): immune GWAS + T-cell sigs
```
Outputs land under `results/real/<dataset>/` as `*.marginal_score.gz` (and, where
`--include_ctrl_score` is used, `*_full.gz`) plus `*.conditional*` tagging scores.

### P2 — DE-overlap simulations
```bash
# 1) generate the simulation inputs (DE h5ad + overlap gene sets):
jupyter nbconvert --to notebook --execute notebooks/02_generate_simulations.ipynb
# 2) score the 5×3×5×2 = 150-job grid (scDRS-FM MAGIC + scDRS none):
bash scripts/09_run_simulations.sh                 # -> results/sim/simulation_data/
```

### P3 — (optional) Null calibration
1,200 null gene sets on TMS_FACS — only needed for the FDP calibration panel:
```bash
bash scripts/10_run_null_calibration.sh            # needs data/gene_sets/gs_files_null/
```

### P4 — Cell-type association tables
Builds the `.scdrs_ct.<biocol>` tables consumed by nb01/06/09/10 (Phase A
re-scores the variants that need control scores; Phase B runs
`scdrs.method.downstream_group_analysis`). Both phases are idempotent:
```bash
bash scripts/13_make_ct_tables.sh                  # or: ... list | rescore-all | downstream-all
```

### P5 — Phenotype gradient decomposition (for nb11)
Runs the package's `run_decompose_gradients.py` on the T-cell / microglia
phenotype scores:
```bash
bash scripts/11_run_phenotype_decomposition.sh     # -> *.marg-marg.decomposition.tsv.gz
```

### P6 — MAGMA gene-level Z matrix (for nb08)
Regenerates the wide gene × trait ZSTAT matrix from the LDSC sumstats using MAGMA
(1000G EUR, NCBI37.3, 10-kb / 10-kb window):
```bash
bash scripts/12_run_magma_geneanalysis.sh          # -> magma_ref/MAGMA_v108_GENE_10_ZSTAT_for_scDRS.txt
```
> The MAGMA binary must live on a filesystem that allows the exec bit (local
> disk). The output filename keeps the paper's `v108` convention; the binary
> actually shipped is MAGMA **v1.10** (see caveats below).

### P7 — Analysis notebooks
The notebooks are already executed with outputs in `notebooks/`. To **re-run** them
after producing the scores above, execute in this order (02 is run in P2):
```
01_problem_fig            05_TMS_Heatmaps           09_soskic_analysis
03_evaluate_simulations   06_watanabe_overlap       10_sea_ad_analysis
                          07_TMSD_TSF_heatmaps      11_phenotype_analysis
                          08_genetic_correlation
```
`04_monocyte_sims.ipynb` is shipped **unrun / optional** (see caveats).

`06_watanabe_overlap.ipynb` reads its inputs entirely through **environment
variables** (no path edits). Set these before executing it:

| Env var | Point at |
|---|---|
| `SCDRSFM_RESULTS_DIR` | `results/ct/tms_facs_magic_ctrl` (has `.scdrs_ct.cell_ontology_class`) |
| `SCDRS_RESULTS_DIR` | `results/ct/tms_facs_none_ctrl` (baseline) |
| `TMS_H5AD_FILE` | `data/subsets_10k/TMS_FACS/TMS_FACS.h5ad` |
| `TMS_CELLTYPE_COLUMN` | `cell_ontology_class` |
| `TMS_CELL_ID_COLUMN` | `cell_id` |
| `TMS_MIN_GENES_PER_CELL` | `250` |

(The notebook also honors `SCDRSFM_RESULTS_DIR_CANDIDATES` /
`SCDRS_RESULTS_DIR_CANDIDATES` colon-lists; the full hook list is in
[`docs/NB_EDIT_MAP.md`](docs/NB_EDIT_MAP.md).)

### P8 — New figures
```bash
python new_figures/make_fig_runtime_memory.py          # runtime + peak-memory profile (Fig A)
python new_figures/make_fig_imputation_none_scatter.py # scDRS-FM vs scDRS concordance (Fig B)
```

---

## 5. Reproduction caveats (read before comparing to the paper)

1. **10k-cell subsampling.** Every real dataset is subsampled to 10,000 cells
   (`seed=0`, `min_genes>=250`) so the whole pipeline runs on one machine.
   Absolute effect sizes and significance counts therefore differ from the
   full-cohort paper numbers, but the qualitative results and method comparisons
   reproduce. Details in [`docs/DATA.md`](docs/DATA.md).
2. **Executed-notebook output paths.** The notebooks are shipped with the outputs
   captured on the reproduction machine, so some **printed** paths in the saved
   cell outputs read `/mnt/shared-workspace/...`. These are recorded stdout from
   the original run and are intentionally **not** rewritten (editing captured
   outputs would misrepresent what actually ran). The notebook *code* is fully
   repo-relative / env-driven — re-executing writes to `results/` under the repo.
3. **MAGMA version label.** `download/download_magma_reference.sh` installs MAGMA
   **v1.10**, but the generated matrix keeps the scDRS naming convention
   `MAGMA_v108_GENE_10_ZSTAT_for_scDRS.txt`. Same file format; only the filename
   label differs.
4. **Genetic correlation (nb08).** The public release ships an `rg` matrix
   (`genet_cor.csv`, 74×74) but **not** the per-pair LDSC logs, so no
   genetic-correlation p-values are available. nb08 loads the provided `rg` values
   and sets p-values to `NaN`, i.e. the `rg` heatmap is shown **without
   significance stars**. The notebook's own unmodified sparse-trait filter then
   drops the all-NaN `PASS_Type_1_Diabetes` row → 73 traits downstream. All logic
   is otherwise unchanged.
5. **scPagwas panels omitted.** scPagwas comparison data is not part of the
   release. nb09 and nb10 guard each scPagwas load with
   `try/except FileNotFoundError` and cleanly skip only those panels; all
   scDRS / scDRS-FM logic and statistics are untouched.
6. **`04_monocyte_sims.ipynb` shipped unrun.** The monocyte simulation panel is
   optional and not required for any headline result; it is included with the
   original outputs but not re-executed in this reproduction.
7. **Simulations use 3 replicates.** The DE-overlap grid uses `reps=0..2` (3
   replicates) rather than the paper's larger replicate count, to keep the
   simulation sweep tractable at 10k scale.

---

## 6. New figures added by this reproduction

* **Fig A — runtime & peak-memory profile** (`new_figures/make_fig_runtime_memory.py`).
  Per-trait and one-time-phase wall-clock timings for the TMS_FACS scDRS-FM run
  (75 traits) plus a real peak-RSS measurement
  (`resource.getrusage(RUSAGE_CHILDREN)`), reproduced from
  `figA_per_trait_timings.csv` + `memprofile_result.json`.
* **Fig B — scDRS-FM vs scDRS concordance** (`new_figures/make_fig_imputation_none_scatter.py`).
  Per-cell and per-trait comparison of MAGIC-imputed scDRS-FM against standard
  scDRS (imputation none) on TMS_FACS, reproduced from the `*.marginal_score.gz`
  files under `results/real/tms_facs/` and `results/ct/tms_facs_none_ctrl/`.

Both scripts are self-contained and regenerate their `.png`/`.svg` + summary CSVs
from the score files.

---

## 7. Reference

This is a reproduction of the scDRS-FM method and paper. The method, its package,
and the underlying data release are the authors'; see `scDRS-FM-main/README.md`
and `scDRS-FM-main/LICENSE`. scDRS-FM builds on **scDRS**
(Zhang et al., *Nature Genetics* 2022) and **MAGIC**
(van Dijk et al., *Cell* 2018).
