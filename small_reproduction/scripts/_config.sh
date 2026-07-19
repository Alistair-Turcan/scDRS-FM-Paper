#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# _config.sh  — shared configuration sourced by every analysis script.
#
# All paths are resolved RELATIVE TO THE REPO ROOT by default, but each can be
# overridden with an environment variable (useful for cluster/scratch layouts).
#
# The original paper ran these as SLURM array jobs calling `run_scdrs+.py`.
# Here they run as plain bash and call the *unmodified* scDRS-FM package entry
# point `run_scdrs_fm.py` (identical CLI). See scripts/run_scdrs_fm.py.
# ---------------------------------------------------------------------------
set -euo pipefail

# Repo root = parent of this scripts/ dir
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export REPO_ROOT="${REPO_ROOT:-$(cd "${SCRIPT_DIR}/.." && pwd)}"

# --- Inputs (override via env as needed) -----------------------------------
export DATA_DIR="${DATA_DIR:-${REPO_ROOT}/data}"                 # from download_figshare_data.sh
export SUBSET_DIR="${SUBSET_DIR:-${DATA_DIR}/subsets_10k}"       # from 00_make_10k_subsets.sh
export GS_SPLIT="${GS_SPLIT:-${DATA_DIR}/gene_sets/gs_split}"    # 75 GWAS gene sets
export TCELL_PHENO="${TCELL_PHENO:-${DATA_DIR}/gene_sets/t_cell_pheno}"
export MICROGLIA="${MICROGLIA:-${DATA_DIR}/gene_sets/microglia}"
export MAGMA_REF="${MAGMA_REF:-${REPO_ROOT}/magma_ref}"

# --- Outputs ---------------------------------------------------------------
export RESULTS="${RESULTS:-${REPO_ROOT}/results}"
export REAL_OUT="${REAL_OUT:-${RESULTS}/real}"
export CT_OUT="${CT_OUT:-${RESULTS}/ct}"
export SIM_OUT="${SIM_OUT:-${RESULTS}/sim/simulation_data}"
export LOG_DIR="${LOG_DIR:-${RESULTS}/logs}"
mkdir -p "${REAL_OUT}" "${CT_OUT}" "${SIM_OUT}" "${LOG_DIR}"

# --- scDRS-FM package entry point (UNMODIFIED) -----------------------------
# Point SCDRS_FM_PY at scDRS-FM-main/run_scdrs_fm.py. Default assumes the
# package sits next to the repo (or vendored under ./scDRS-FM-main for convenience).
export SCDRS_FM_HOME="${SCDRS_FM_HOME:-${REPO_ROOT}/scDRS-FM-main}"
export SCDRS_FM_PY="${SCDRS_FM_PY:-${SCDRS_FM_HOME}/run_scdrs_fm.py}"
export PYTHON="${PYTHON:-python}"

# Make scdrs_fm importable when calling the script directly
export PYTHONPATH="${SCDRS_FM_HOME}:${PYTHONPATH:-}"

# All datasets are subsampled to 10k cells before scoring (fixed seed=0).
# See 00_make_10k_subsets.sh / docs/DATA.md.

run_scdrs_fm () {
  echo "[run_scdrs_fm] ${PYTHON} ${SCDRS_FM_PY} $*"
  "${PYTHON}" "${SCDRS_FM_PY}" "$@"
}
