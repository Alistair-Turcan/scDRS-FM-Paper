#!/usr/bin/env python3
"""
Gradient-decomposition pipeline step for notebook 11 (phenotype analysis).

nb11 READS pre-computed `{trait}.marg-marg.decomposition.tsv.gz` files but does NOT
produce them. They come from the scDRS-FM package's decomposition stage
(`scdrs_fm.decompose_gradients.run_decomposition`, CLI: run_decompose_gradients.py),
which the paper ran ad-hoc (no slurm in the repo). This driver runs that exact
package function -- the package is NOT modified -- with the user-approved config:

  * Phenotype predictors  = the 24 FUNCTIONAL signatures curated in nb09
                            (microglia dataset uses its 5 microglia signatures).
  * Trait/outcome scores  = MAGIC-imputed scores.
  * Cells                 = ALL cells in each 10k subset (no CD4 subsetting;
                            `cd4_*` in the nb11 figure names is a label only).

Outputs are written into each dataset's pheno_dir (== the dir nb11 reads:
RESULTS/real/{soskic_tcell,nathan_tcell,canogamez_tcell,sea_ad_micro}).

Deterministic (n_controls=1000, the package's MC seed handling). Idempotent:
skips a dataset if all its expected decomposition files already exist.
"""
from __future__ import annotations
import os
import shutil
import sys
import time
from pathlib import Path

BASE = Path(os.environ.get("REPO_ROOT", "."))
sys.path.insert(0, os.environ.get("SCDRS_FM_HOME", str(BASE / "scDRS-FM-main")))
from scdrs_fm.decompose_gradients import run_decomposition  # noqa: E402

DATA = Path(os.environ.get("SUBSET_DIR", str(BASE / "data" / "subsets_10k")))
RESULTS = Path(os.environ.get("RESULTS", str(BASE / "results")))

# 24 functional phenotype signatures (nb09 PHENOTYPE_TRAITS) -- user decision.
SIG24 = [
    "Metallothionein", "Translation", "IL10-IL19", "OX40-EBI3", "CD172a-MERTK",
    "TIMD4-TIM3", "BCL2-FAM13A", "IEG", "SOX4-TOX2", "NME1-FABP5", "IEG3",
    "RGCC-MYADM", "Exhaustion", "ISG", "Cytotoxic", "CD40LG-TXNIP", "Mito",
    "HLA", "Heatshock", "IEG2", "Cytoskeleton", "CTLA4-CD38", "Multi-Cytokine",
    "ICOS-CD38",
]

# 10 immune GWAS traits (nb11 PHENOTYPE_TRAITS for the T-cell datasets).
IMMUNE_TRAITS = [
    "PASS_CD_deLange2017", "PASS_UC_deLange2017", "PASS_IBD_deLange2017",
    "PASS_Celiac", "PASS_Rheumatoid_Arthritis", "UKB_460K.disease_AID_ALL",
    "UKB_460K.disease_ASTHMA_DIAGNOSED", "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED",
    "UKB_460K.disease_RESPIRATORY_ENT", "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP",
]

# 5 microglia signatures + 2 microglia GWAS traits (nb11 MICROGLIA_DATASET).
MICRO_SIGS = ["CRM_gs", "DAM_gs", "HLA_gs", "HM_gs", "IRM_gs"]
MICRO_TRAITS = ["PASS_Parkinsons23andMe_Corces2020", "PASS_Alzheimers_Jansen2019"]

# (name, adata, trait_dir, pheno_dir==out_dir, traits, phenotypes)
CONFIGS = [
    ("soskic", DATA / "Soskic" / "soskic_100k.h5ad",
     RESULTS / "ct" / "soskic_immune_magic_ctrl", RESULTS / "real" / "soskic_tcell",
     IMMUNE_TRAITS, SIG24),
    ("nathan", DATA / "Nathan" / "raw.h5ad",
     RESULTS / "real" / "nathan_immune", RESULTS / "real" / "nathan_tcell",
     IMMUNE_TRAITS, SIG24),
    ("canogamez", DATA / "Cano_Gamez" / "obj_raw.h5ad",
     RESULTS / "real" / "canogamez_immune", RESULTS / "real" / "canogamez_tcell",
     IMMUNE_TRAITS, SIG24),
    ("microglia", DATA / "SEA_AD" / "combined_healthy_filtered.h5ad",
     RESULTS / "ct" / "sea_ad_brain_magic_ctrl", RESULTS / "real" / "sea_ad_micro",
     MICRO_TRAITS, MICRO_SIGS),
]


def already_done(out_dir: Path, traits) -> bool:
    return all((out_dir / f"{t}.marg-marg.decomposition.tsv.gz").exists() for t in traits)


def stage_trait_dir(name: str, trait_dir: Path, traits) -> Path:
    """The decomposition's outcome/trait loaders require 1000 `ctrl_norm_score_*`
    columns in both the marginal and conditional trait tables. The thinned
    `{trait}.marginal_score.gz` / `{trait}.conditional.tagging_score.gz` files do
    NOT carry ctrl columns, but the `*_full.gz` variants produced by the same
    scoring run DO (verified: 1000 ctrl cols, correct schema). Rather than
    re-score or mutate the canonical thinned files (read by other notebooks),
    stage a local dir where the `_full` files are renamed to the names
    `run_decomposition` expects. No package modification; no re-scoring.
    """
    stage = Path(os.environ.get("LOCAL_TMP", "/tmp")) / "decomp_stage" / name
    stage.mkdir(parents=True, exist_ok=True)
    for t in traits:
        pairs = [
            (trait_dir / f"{t}.marginal_score_full.gz", stage / f"{t}.marginal_score.gz"),
            (trait_dir / f"{t}.conditional_score_full.gz", stage / f"{t}.conditional.tagging_score.gz"),
        ]
        for src, dst in pairs:
            if not src.exists():
                raise FileNotFoundError(f"[{name}] required full-ctrl file missing: {src}")
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            # symlink to avoid copying multi-MB ctrl tables; pandas reads through it
            os.symlink(src, dst)
    return stage


def main() -> None:
    only = set(sys.argv[1:])  # optional: restrict to named datasets
    for name, adata, trait_dir, pheno_dir, traits, phenos in CONFIGS:
        if only and name not in only:
            continue
        print(f"\n{'='*70}\n[{name}] decomposition\n{'='*70}", flush=True)
        if already_done(pheno_dir, traits):
            print(f"[{name}] all {len(traits)} decomposition files already exist -- skip", flush=True)
            continue
        t0 = time.perf_counter()
        staged = stage_trait_dir(name, trait_dir, traits)
        print(f"[{name}] staged {len(traits)} full-ctrl trait tables -> {staged}", flush=True)
        run_decomposition(
            adata_file=str(adata),
            trait_dir=str(staged),
            pheno_dir=str(pheno_dir),
            out_dir=str(pheno_dir),           # nb11 reads decomposition from pheno_dir
            traits=list(traits),
            phenotypes=list(phenos),
            additional_phenotypes=[],
            cell_type_col="",                 # all cells (user decision)
            cell_type="",
            n_controls=1000,
        )
        # verify
        made = [t for t in traits if (pheno_dir / f"{t}.marg-marg.decomposition.tsv.gz").exists()]
        print(f"[{name}] wrote {len(made)}/{len(traits)} decomposition files "
              f"in {time.perf_counter()-t0:,.1f}s", flush=True)
        if len(made) != len(traits):
            missing = [t for t in traits if t not in made]
            print(f"[{name}] WARNING missing: {missing}", flush=True)


if __name__ == "__main__":
    main()
