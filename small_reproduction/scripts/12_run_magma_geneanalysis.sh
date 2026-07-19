#!/usr/bin/env bash
# 12_run_magma_geneanalysis.sh  — regenerate the gene-level Z matrix for nb08.
#
# Pipeline (MAGMA >= v1.08 semantics; 1000G Phase 3 EUR; NCBI37.3; window 10,10 kb):
#   1. sumstats/*.sumstats.gz  --(Z -> two-sided P)-->  magma_ref/pval/<trait>.{pval,N}
#   2. magma --bfile g1000_eur --pval <trait>.pval N=<N> --gene-annot ...w10.genes.annot
#        -> magma_ref/out/<trait>.genes.out            (per trait)
#   3. assemble_zstat.py  -> MAGMA_v108_GENE_10_ZSTAT_for_scDRS.txt (wide gene x trait ZSTAT)
#
# Prereq: download/download_magma_reference.sh (binary + 1000G EUR + gene.loc + annot).
# The MAGMA binary must run from local disk (exec bit); copy magma_ref to local
# scratch if magma_ref sits on a network mount.
#
# Usage:  bash scripts/12_run_magma_geneanalysis.sh
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"

SUMSTATS="${SUMSTATS:-${DATA_DIR}/sumstats}"
PVAL_DIR="${MAGMA_REF}/pval"
OUT_DIR="${MAGMA_REF}/out"
BIN="${MAGMA_BIN:-${MAGMA_REF}/magma}"
BFILE="${MAGMA_BFILE:-${MAGMA_REF}/g1000_eur/g1000_eur}"
ANNOT="${MAGMA_ANNOT:-${MAGMA_REF}/g1000_eur_NCBI37_w10.genes.annot}"
mkdir -p "${PVAL_DIR}" "${OUT_DIR}"

echo "[12/1] sumstats -> pval/N ..."
"${PYTHON}" "$(dirname "${BASH_SOURCE[0]}")/prep_magma_pval.py" "${SUMSTATS}" "${PVAL_DIR}"

echo "[12/2] MAGMA gene analysis per trait ..."
if [[ ! -x "${BIN}" ]]; then
  echo "ERROR: MAGMA binary not executable at ${BIN}. Run download_magma_reference.sh,"
  echo "       and ensure magma_ref is on a filesystem that allows the exec bit (local disk)."
  exit 2
fi
shopt -s nullglob
for pv in "${PVAL_DIR}"/*.pval; do
  t="$(basename "${pv}" .pval)"
  if [[ -f "${OUT_DIR}/${t}.genes.out" ]]; then echo "  [skip] ${t}"; continue; fi
  N="$(cat "${PVAL_DIR}/${t}.N")"
  echo "  [magma] ${t}  N=${N}"
  "${BIN}" --bfile "${BFILE}" --pval "${pv}" N="${N}" \
     --gene-annot "${ANNOT}" --out "${OUT_DIR}/${t}" > "${OUT_DIR}/${t}.log" 2>&1 || \
     echo "  [FAIL] ${t} (see ${OUT_DIR}/${t}.log)"
done

echo "[12/3] Assembling wide ZSTAT matrix ..."
export MAGMA_REF MAGMA_OUT="${OUT_DIR}"
"${PYTHON}" "$(dirname "${BASH_SOURCE[0]}")/assemble_zstat.py"
echo "[12] Done -> ${MAGMA_REF}/MAGMA_v108_GENE_10_ZSTAT_for_scDRS.txt"
