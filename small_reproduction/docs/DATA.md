# Data & the 10k-cell subsampling design

This document describes the datasets used by the reproduction and the single
substantive deviation from the paper: subsampling every real dataset to 10,000
cells. All of this lives in repo-tooling only — the scDRS-FM package under
`scDRS-FM-main/` is untouched.

## 1. Source data (figshare article 33000602)

`download/download_figshare_data.sh` pulls the scDRS-FM data release:

| Dataset | h5ad | cov | species | full cells |
|---|---|---|---|---|
| TMS_FACS | `TMS_FACS.h5ad` | `TMS_FACS.cov` | mouse | 110,824 |
| TS_FACS | `ts_facs.h5ad` | `ts.cov` | human | 483,152 |
| TMS_Droplet | `TMS_Droplet.h5ad` | `tms_droplet.cov` | mouse | 245,389 |
| SEA_AD | `combined_healthy_filtered.h5ad` | `combined_healthy_filtered.cov` | human | 168,478 |
| Soskic | `soskic_100k.h5ad` | `soskic.cov` | human | 651,650 |
| Nathan | `raw.h5ad` | `raw.cov` | human | 500,089 |
| Cano_Gamez | `obj_raw.h5ad` | `canogamez.cov` | human | 43,112 |
| Braun | `human_dev_layers_100k.h5ad` | `human_dev_layers_100k.cov` | human | 100,000 |

Plus gene sets and annotations:

* `GS_files_75.zip` → `gs_split/` — 75 MAGMA GWAS gene sets (one file per trait).
* `Microglia_phenotype_genesets.zip` → `microglia/` — `HM_gs`, `DAM_gs`,
  `CRM_gs`, `IRM_gs`, `HLA_gs` (5 signatures).
* `T_cell_phenotype_genesets.zip` → `t_cell_pheno/` — 52 T-cell phenotype signatures.
* `genet_cor.csv` — LDSC genetic-correlation (`rg`) matrix (74×74). **No** per-pair
  p-values are released (see nb08 caveat in the top-level README).
* `*.sumstats.gz` → `sumstats/` — 75 LDSC-format GWAS summary statistics
  (`SNP A1 A2 Z N`), consumed by the MAGMA step.

The full inventory is in [`../download/manifest.tsv`](../download/manifest.tsv).

## 2. Why subsample to 10k cells

The paper scored full cohorts (tens to hundreds of thousands of cells per
dataset) on a SLURM cluster. To make the **entire** pipeline reproducible
end-to-end on a single workstation — including MAGIC imputation, which is the
memory/compute bottleneck of scDRS-FM — every real dataset is subsampled to
exactly **10,000 cells** before scoring.

This is the only substantive deviation from the paper. It is implemented purely
in `scripts/make_10k_subsets.py` (repo tooling); the scDRS-FM package and the CLI
flags used for each dataset are identical to the originals.

## 3. Exactly how the subset is drawn (`scripts/00_make_10k_subsets.sh`)

For each dataset:

1. **Eligibility filter.** The scDRS-FM pipeline (`run_scdrs_fm.py --flag_filter`)
   applies `sc.pp.filter_cells(min_genes=250)` *before* scoring. To guarantee that
   a full 10,000 cells actually reach MAGIC + scoring (rather than fewer after
   in-pipeline filtering), we subsample from the pool of cells that **already pass
   `n_genes >= 250`**. `n_genes` is read from the covariate file when available and
   otherwise computed from `X` in memory-bounded chunks.
2. **Random draw.** From the eligible pool, draw exactly `TARGET_N = 10,000` cells
   without replacement using `numpy.random.default_rng(seed=0)`
   (**fixed `seed=0`** for full determinism). If fewer than 10,000 cells are
   eligible, all eligible cells are kept.
3. **Preserve raw data.** The raw `X` matrix and the **full gene set** are
   preserved (gene-level filtering still happens in-pipeline on the 10k subset,
   exactly as in the original per-dataset run). The covariate file is reindexed to
   the selected barcodes and written alongside.
4. **Braun name reconciliation.** figshare ships
   `human_dev_layers_100k.{h5ad,cov}`, but some original SLURM scripts reference
   `*_gene_symbols.h5ad` / `*_processed.cov`. The subset step writes **both** names
   (copies) so either path resolves.

Outputs land in `data/subsets_10k/<Dataset>/` together with a
`subset_summary.json` recording `seed`, `target_n`, `min_genes`, and per-dataset
cell counts (`n_total`, `n_eligible`, `n_selected`, `n_genes_full`).

## 4. Per-dataset scoring flags (unchanged from the paper)

The subsampling does not change any scoring flag. Each `scripts/0X_run_*.sh`
reproduces the exact per-dataset CLI from the original SLURM jobs:

| Dataset | species | filter | raw count | imputation | ctrl scores | trait set(s) |
|---|---|---|---|---|---|---|
| TMS_FACS | mouse | yes | yes | **magic** | no | 75 GWAS |
| TS_FACS | human | yes | yes | none | yes | 75 GWAS |
| TMS_Droplet | mouse | yes | yes | none | yes | 75 GWAS |
| SEA_AD | human | yes | no | **magic** | (ct step) | 8 brain GWAS + microglia |
| Soskic | human | yes | yes | **magic** | (ct step) | 15 immune GWAS + 52 T-cell sigs |
| Nathan | human | yes | yes | **magic** | yes | immune GWAS + T-cell sigs |
| Cano_Gamez | human | yes | yes | **magic** | yes | immune GWAS + T-cell sigs |
| Braun | human | yes | yes | **magic** | yes | 8 brain GWAS + 5 microglia |
| Simulations | human | yes | yes | magic (scDRS-FM) / none (scDRS) | — | DE-overlap grid |
| Null calib. | human | yes | yes | **magic** | — | 1,200 null gene sets |

"magic" = scDRS-FM (MAGIC imputation, `random_state=0`); "none" = standard scDRS.
The `.scdrs_ct.<biocol>` cell-type tables consumed by nb01/06/09/10 are produced
by `scripts/13_make_ct_tables.sh` (which re-scores with `--include_ctrl_score`
where needed and runs `scdrs.method.downstream_group_analysis`).

## 5. Consequence for the results

Because scores are computed on 10k-cell subsets, **absolute** effect sizes and
counts of FDR-significant cells differ from the full-cohort paper numbers. The
**qualitative** findings and the method comparisons (scDRS-FM vs scDRS) reproduce.
Every notebook is shipped executed so the reproduced numbers are visible without
rerunning; see the executed-output-path caveat in the top-level README.
