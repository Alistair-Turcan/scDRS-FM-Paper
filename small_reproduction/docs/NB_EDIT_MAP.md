# Notebook path-edit map (P5) — scDRS-FM reproduction @ 10k scale

All notebooks executed from a working dir with these output/result trees available:
- REAL scores:  /mnt/shared-workspace/scdrsfm/results/real/<job_id>/
- CT rescored:  /mnt/shared-workspace/scdrsfm/results/ct/<job_id>/  (has .scdrs_ct.<biocol>)
- SIM scores:   /mnt/shared-workspace/scdrsfm/results/sim/simulation_data/  (= nb03 BASE_DIR)
- Subsets:      /mnt/shared-workspace/scdrsfm/data/subsets_10k/<DS>/
- Gene sets:    /mnt/shared-workspace/scdrsfm/data/gene_sets/{gs_split,t_cell_pheno,microglia}
- MAGMA zstat:  /mnt/shared-workspace/scdrsfm/magma_ref/MAGMA_v108_GENE_10_ZSTAT_for_scDRS.txt (regenerated)
- genet_cor:    /mnt/shared-workspace/scdrsfm/data/genet_cor.csv (rg only, 74x74; NO pval matrix available)

INDEP_CELLS: RESOLVED — nb05/07/09/10 each contain a function that READS marginal+conditional
scores and WRITES indep_cells/<...>.gz. None REQUIRE pre-existing indep_cells. Leave INDEP_CELLS_DIR
pointing at a writable local path (default relative path is fine; it mkdirs).

## nb01_problem_fig  — scDRS baseline discovered-CT boxplot on TMS_FACS
- RESULTS_DIR = Path("tms_data/scdrs_results")  -> results/ct/tms_facs_none_ctrl
  (suffix .scdrs_ct.cell_ontology_class ; baseline = imputation NONE)
- Reads .scdrs_ct.cell_ontology_class for BRAIN/IMMUNE trait lists (defined inline).

## nb02_generate_simulations  — DONE already (produced gs + DE h5ad). No exec needed at P5.

## nb03_evaluate_simulations  — single BASE_DIR edit
- BASE_DIR -> results/sim/simulation_data (contains scdrs+_results, scdrs_results,
  *_cell_pct, *_causal_genes, scdrs+_results_knn + gs_de_overlap/gs_cell_pct/gs_causal_genes
  + adata/ + predictions*.csv). Run wire_nb03_inputs.py first to stage the gs/predictions/adata.

## nb04_monocyte_sims  — SHIP UNRUN (user decision #2). No exec.

## nb05_TMS_Heatmaps  — scDRS-FM TMS_FACS heatmaps + writes indep_cells
- H5AD_FILE   = "../scdrs_simpleby/tms_data/.../TMS_FACS.h5ad" -> subsets_10k/TMS_FACS/TMS_FACS.h5ad
- GS_DIR      = Path("../scdrs_simpleby/gs_split")             -> data/gene_sets/gs_split
- RESULTS_DIR = Path("../scdrs_simpleby/august_all/scdrs+_results_gram3") -> results/real/tms_facs
  (tms_facs = magic no-ctrl; nb05 only needs marginal+conditional, which tms_facs has)
- INDEP_CELLS_DIR = Path("indep_cells/tms_facs")  -> leave (writable, notebook writes it)

## nb06_watanabe_overlap  — ENV-VAR HOOKS (no code edit). scDRS-FM vs scDRS comparison
- Set env vars before running:
  SCDRSFM_RESULTS_DIR = results/ct/tms_facs_magic_ctrl   (needs .scdrs_ct.cell_ontology_class -> FM)
  SCDRS_RESULTS_DIR   = results/ct/tms_facs_none_ctrl    (baseline .scdrs_ct.cell_ontology_class)
  TMS_H5AD_FILE       = subsets_10k/TMS_FACS/TMS_FACS.h5ad
  (12 hooks total; confirm full list from nb06 source at exec time)

## nb07_TMSD_TSF_heatmaps  — TMS_Droplet + TS_FACS scDRS-FM heatmaps; writes indep_cells
- h5ad_file=Path("tms_droplet/TMS_Droplet.h5ad") -> subsets_10k/TMS_Droplet/TMS_Droplet.h5ad
- results_dir=Path("ct_validations/tms_droplet/scdrs+_results") -> results/real/tms_droplet
- h5ad_file=Path("ts_data/ts_facs.h5ad") -> subsets_10k/TS_FACS/ts_facs.h5ad
- results_dir=Path("ct_validations/ts_facs/scdrs+_results") -> results/real/ts_facs
- STRICT_INPUT_FILES = True (keep; all inputs exist)
- NOTE: tms_droplet + ts_facs were scored with imputation=none + ctrl (per real_jobs). nb07 reads
  marginal+conditional.tagging_score.gz which both have.

## nb08_genetic_correlation  — MAGMA zstat + genet_cor wiring + scDRS-FM scores
- MAGMA_GENE_ZSTAT_FILE -> magma_ref/MAGMA_v108_GENE_10_ZSTAT_for_scDRS.txt
- GS_SPLIT_DIR = Path("../scdrs_simpleby/gs_split") -> data/gene_sets/gs_split
- Score dirs:
    DatasetConfig("tms_facs", "august_all/scdrs+_results_gram3") -> results/real/tms_facs
    DatasetConfig("ts_facs",  "ct_validations/ts_facs/scdrs+_results") -> results/real/ts_facs
    DatasetConfig("tms_droplet","ct_validations/tms_droplet/scdrs+_results") -> results/real/tms_droplet
- GENET_COR WIRING (decision #10 + USER DECISION 2026-07-18 = KEEP RG, DROP STARS):
  In cell 24, REPLACE these two lines:
      genet_cor, genet_cor_pval = build_genetic_correlation_matrices(ALL_TRAITS, LDSC_COR_DIR)
      genet_cor, genet_cor_pval = drop_sparse_genetic_traits(genet_cor, genet_cor_pval, max_missing=5)
  WITH (load provided rg csv; pval all-NaN so no FDR stars; KEEP drop_sparse unchanged):
      # LDSC per-pair logs are not in the public release; load the provided rg matrix instead
      # and set p-values to NaN (no genetic significance stars). Everything else unchanged.
      _rg_raw = pd.read_csv(GENET_COR_CSV, index_col=0)
      genet_cor = _rg_raw.reindex(index=ALL_TRAITS, columns=ALL_TRAITS)
      genet_cor_pval = pd.DataFrame(np.nan, index=ALL_TRAITS, columns=ALL_TRAITS)
      genet_cor, genet_cor_pval = drop_sparse_genetic_traits(genet_cor, genet_cor_pval, max_missing=5)
  Also add to the Input-paths cell (cell 3):  GENET_COR_CSV = Path("<abs>/data/genet_cor.csv")
  VALIDATED OUTCOME: provided rg csv has 74 traits INCLUDING Parkinsons but MISSING PASS_Type_1_Diabetes;
  ALL_TRAITS (gs_split minus TRAITS_TO_DROP=Parkinsons) has 74 INCLUDING T1D. After reindex+drop_sparse,
  T1D (all-NaN row) is dropped by the notebook's OWN unmodified sparse filter -> 73 traits downstream.
  All 16 SUBSET_TRAITS still present (T1D not among them). rg heatmap colors fully present; no stars.
  build_gene_level_correlation(MAGMA, TRAITS) uses the 73-trait set; MAGMA matrix (74, T1D present,
  Parkinsons absent) intersects fine. TRAITS_TO_DROP, drop_sparse_genetic_traits, plot logic UNCHANGED.

## nb09_soskic_analysis  — t-cell; scDRS-FM vs scDRS vs scPagwas; reads .scdrs_ct
- Uses first_existing_path candidate lists. Point candidates (or place files) at:
    SCDRSFM_SOSKIC_RESULTS  <- results/ct/soskic_immune_magic_ctrl  (.marginal_score.gz.scdrs_ct.Cell_population)
    SCDRS_SOSKIC_RESULTS    <- results/ct/soskic_immune_none_ctrl   (baseline)
    PHENO_SOSKIC_RESULTS    <- results/real/soskic_tcell            (t_cell_pheno marginal scores)
    T_CELL_PHENO_DIR        <- data/gene_sets/t_cell_pheno
    soskic h5ad             <- subsets_10k/Soskic/soskic_100k.h5ad
    Nathan  scdrsfm         <- results/ct/nathan_*  (immune ct done; tcell magic ctrl for pheno)
    Cano_Gamez scdrsfm      <- results/ct/canogamez_*
- SCPAGWAS_SOSKIC_DIR -> GUARD (decision #1): skip scPagwas panels cleanly.
- Writes indep_cells (function present). No pre-existing needed.

## nb10_sea_ad_analysis  — SEA_AD + Braun; scDRS-FM vs scDRS vs scPagwas; reads .scdrs_ct
- load_and_preprocess_adata(path, normalize=False): only uses obs labels + counts (NOT expression).
  Raw-count subsets are fine.
- SEA_AD_H5AD            <- subsets_10k/SEA_AD/combined_healthy_filtered.h5ad
- SEA_AD_SCDRSFM_DIR     <- results/ct/sea_ad_brain_magic_ctrl (.scdrs_ct.Subclass)
- SEA_AD_SCDRS_DIR       <- results/ct/sea_ad_brain_none_ctrl (baseline)
- SEA_AD_PHENOTYPE_SCORE_DIR <- results/real/sea_ad_micro (HM_gs.marginal_score.gz etc.)
- BRAUN_H5AD             <- subsets_10k/Braun/human_dev_layers_100k.h5ad (has CellClass+Region obs)
- BRAUN_SCDRSFM_DIR      <- results/ct/braun_micro  OR results/real/braun_micro (both have full+ct)
- SEA_AD_SCPAGWAS_DIR / SEA_AD_SCPAGWAS_AD_DIR -> GUARD (decision #1): skip scPagwas panels cleanly.
- Writes indep_cells. No pre-existing needed.

## nb11_phenotype_analysis  — CD4 phenotype marginal scores (Soskic/Cano/Nathan/SEA_AD micro)
- results_dir=Path("ct_validations/soskic/scdrs+_results")   -> results/real/soskic_tcell
- results_dir=Path("ct_validations/canogamez/scdrs+_results")-> results/real/canogamez_tcell
- results_dir=Path("ct_validations/nathan/scdrs+_results")   -> results/real/nathan_tcell
- results_dir=Path("sea_ad_final/healthy_results_filtered_martin") -> results/real/sea_ad_micro
- Reads .marginal_score.gz directly (t_cell_pheno gene set scores). No .scdrs_ct needed.

## scPagwas guard strategy (nb09, nb10)
scPagwas data dirs never exist here. Guard each scPagwas LOAD/BUILD call so a missing dir
skips that panel (set the scpagwas table to None / empty and skip the corresponding plot section),
WITHOUT touching scDRS/scDRS-FM logic or plot stats. Minimal try/except FileNotFoundError -> skip.
