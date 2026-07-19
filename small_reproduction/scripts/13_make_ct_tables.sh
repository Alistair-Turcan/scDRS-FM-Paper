#!/usr/bin/env bash
# 13_make_ct_tables.sh  — build the .scdrs_ct.<biocol> cell-type association
# tables consumed by notebooks 01, 06, 09 and 10.
#
# What it does (see scripts/ct_driver.py for the full matrix):
#   Phase A "rescore"  : re-scores the datasets/variants that were NOT scored
#                        with --include_ctrl_score, producing *.marginal_score_full.gz
#                        under results/ct/<id>/. (magic = scDRS-FM, none = std scDRS.)
#   Phase B "downstream": runs scdrs.method.downstream_group_analysis on every
#                        *.marginal_score_full.gz and writes the per-cell-type
#                        association tables (suffix .scdrs_ct.<biocol> or
#                        .marginal_score.gz.scdrs_ct.<biocol> for nb09).
#
# The 3 datasets already scored with control (nathan_immune, canogamez_immune,
# braun_micro) skip Phase A and read their full scores from results/real/<id>/.
#
# Both phases are idempotent (per-job .DONE / .DOWNSTREAM_DONE markers), so the
# script can be safely re-run; finished jobs are skipped.
#
# Usage:
#   bash scripts/13_make_ct_tables.sh                 # rescore-all then downstream-all
#   bash scripts/13_make_ct_tables.sh list            # show job status
#   bash scripts/13_make_ct_tables.sh rescore-all
#   bash scripts/13_make_ct_tables.sh downstream-all
#   bash scripts/13_make_ct_tables.sh rescore <id>...
#   bash scripts/13_make_ct_tables.sh downstream <result_dir>...
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"
export REPO_ROOT DATA_DIR SUBSET_DIR RESULTS CT_OUT LOG_DIR SCDRS_FM_HOME
export LOCAL_TMP="${LOCAL_TMP:-/tmp/scdrsfm}"

CT="$(dirname "${BASH_SOURCE[0]}")/ct_driver.py"

if [[ $# -eq 0 ]]; then
  echo "[13] Phase A: rescore-all (control-augmented scoring) ..."
  "${PYTHON}" "${CT}" rescore-all
  echo "[13] Phase B: downstream-all (cell-type association tables) ..."
  "${PYTHON}" "${CT}" downstream-all
  echo "[13] Done. Cell-type tables written under ${CT_OUT} and results/real/*."
else
  "${PYTHON}" "${CT}" "$@"
fi
