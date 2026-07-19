#!/usr/bin/env bash
# 09_run_simulations.sh  — DE-overlap simulation grid.
# Original: slurms/118_sims.slurm (a 0-149 SLURM array). Reproduced here as a
# serial/loop bash job over the identical grid:
#
#   clusters = 0..4  (5)      reps = 0..2  (3)      src_n = 1..5  (5)
#   methods  = { magic (scDRS-FM) , none (scDRS) }  (2)     overlap = 50
#   => 5 * 3 * 5 * 2 = 150 scoring jobs
#
# Inputs come from notebook 02 (02_generate_simulations.ipynb), which writes:
#   data/simulation_data/adata/TMS_FACS_10k_DE.h5ad
#   data/simulation_data/gs_de_overlap/TMS_FACS_{c}_{r}_50_src{s}.gs
# Outputs (consumed by nb03) go under:
#   results/sim/simulation_data/scdrs+_results/{c}/{r}/src{s}   (scDRS-FM, MAGIC)
#   results/sim/simulation_data/scdrs_results/{c}/{r}/src{s}    (scDRS, none)
#
# NOTE: run notebooks/02_generate_simulations.ipynb FIRST to create the inputs.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"

SIM_BASE="${SIM_BASE:-${DATA_DIR}/simulation_data}"     # inputs (from nb02)
GS_DE="${SIM_BASE}/gs_de_overlap"
H5AD="${SIM_BASE}/adata/TMS_FACS_10k_DE.h5ad"
COV="${SUBSET_DIR}/TMS_FACS/TMS_FACS.cov"
OVERLAP=50

n_cluster=5; n_rep=3; n_src=5
echo "[09] DE-overlap grid: ${n_cluster} clusters x ${n_rep} reps x ${n_src} src x 2 methods = $((n_cluster*n_rep*n_src*2)) jobs"

for cluster in $(seq 0 $((n_cluster-1))); do
  for rep in $(seq 0 $((n_rep-1))); do
    for src in $(seq 1 ${n_src}); do
      trait="TMS_FACS_${cluster}_${rep}_${OVERLAP}_src${src}.gs"
      # (0) scDRS-FM (MAGIC)
      out_fm="${SIM_OUT}/scdrs+_results/${cluster}/${rep}/src${src}"; mkdir -p "${out_fm}"
      run_scdrs_fm "${H5AD}" "${COV}" "${out_fm}" "${GS_DE}" "${trait}" \
        --h5ad_species human --flag_filter --flag_raw_count --imputation magic
      # (1) scDRS (none)
      out_std="${SIM_OUT}/scdrs_results/${cluster}/${rep}/src${src}"; mkdir -p "${out_std}"
      run_scdrs_fm "${H5AD}" "${COV}" "${out_std}" "${GS_DE}" "${trait}" \
        --h5ad_species human --flag_filter --flag_raw_count --imputation none
    done
  done
done
echo "[09] Done -> ${SIM_OUT}/{scdrs+_results,scdrs_results}"
echo "[09] Next: run notebooks/03_evaluate_simulations.ipynb"
