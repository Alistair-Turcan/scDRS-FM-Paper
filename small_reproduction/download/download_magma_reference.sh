#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# download_magma_reference.sh
#
# Fetches the MAGMA binary + reference panel needed to regenerate the
# gene-level Z-statistic matrix used by notebook 08 (genetic correlation).
#
#   - MAGMA v1.10 (linux, statically linked)   from https://cncr.nl/research/magma/
#   - NCBI37.3 gene locations (GRCh37, Entrez)
#   - 1000 Genomes Phase 3 European reference (build 37): g1000_eur.{bed,bim,fam,synonyms}
#
# The output ZSTAT matrix is named MAGMA_v108_GENE_10_ZSTAT_for_scDRS.txt to
# match the scDRS filename convention (window 10kb). MAGMA >= v1.08 semantics.
#
# USAGE:
#   bash download/download_magma_reference.sh [DEST_DIR]
# DEST_DIR defaults to ./magma_ref
#
# NOTE: CNCR/MAGMA share links occasionally rotate. If a URL 404s, get the
# current links from https://cncr.nl/research/magma/ and update below.
# The MAGMA binary MUST be executed from a local (non-network) filesystem so
# the exec bit is honored.
# ---------------------------------------------------------------------------
set -euo pipefail

DEST="${1:-./magma_ref}"
mkdir -p "${DEST}/bin" "${DEST}/gene_loc" "${DEST}/g1000_eur"

MAGMA_URL="https://vu.data.surfsara.nl/index.php/s/lxDgt2dNdNr6DYt/download"          # magma_v1.10 linux static
GENELOC_URL="https://vu.data.surfsara.nl/index.php/s/Pj2orwuF2JYyKxq/download"        # NCBI37.3.gene.loc
G1000_URL="https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download"          # g1000_eur (build 37)

fetch () {  # url  outfile
  if [[ -f "$2" ]]; then echo "  [skip] $(basename "$2")"; return; fi
  echo "  [get ] $(basename "$2")"
  curl -sSL -C - -o "$2" "$1"
}

echo "[magma] Downloading MAGMA binary + reference..."
fetch "${MAGMA_URL}"   "${DEST}/bin/magma_v1.10.zip"
fetch "${GENELOC_URL}" "${DEST}/gene_loc/NCBI37.3.gene.loc.zip"
fetch "${G1000_URL}"   "${DEST}/g1000_eur/g1000_eur.zip"

echo "[magma] Unzipping..."
( cd "${DEST}/bin"       && unzip -n -q magma_v1.10.zip           && chmod +x magma* 2>/dev/null || true )
( cd "${DEST}/gene_loc"  && unzip -n -q NCBI37.3.gene.loc.zip     2>/dev/null || true )
( cd "${DEST}/g1000_eur" && unzip -n -q g1000_eur.zip             2>/dev/null || true )

# Locate the static binary and expose it as ${DEST}/magma
BIN="$(find "${DEST}/bin" -maxdepth 2 -type f -name 'magma*' ! -name '*.zip' | head -1 || true)"
if [[ -n "${BIN}" ]]; then
  cp -f "${BIN}" "${DEST}/magma"
  chmod +x "${DEST}/magma" || true
  echo "[magma] Binary staged at ${DEST}/magma"
  "${DEST}/magma" --version || true
fi

echo "[magma] Precomputing gene annotation (window 10,10 kb)..."
BIM="$(find "${DEST}/g1000_eur" -name 'g1000_eur.bim' | head -1 || true)"
GLOC="$(find "${DEST}/gene_loc" -name 'NCBI37.3.gene.loc' | head -1 || true)"
if [[ -n "${BIM}" && -n "${GLOC}" && -x "${DEST}/magma" ]]; then
  PREFIX="${BIM%.bim}"
  "${DEST}/magma" --annotate window=10,10 \
      --snp-loc "${BIM}" \
      --gene-loc "${GLOC}" \
      --out "${DEST}/g1000_eur_NCBI37_w10" || true
  echo "[magma] Annotation -> ${DEST}/g1000_eur_NCBI37_w10.genes.annot"
fi

echo "[magma] Done. Next: bash scripts/12_run_magma_geneanalysis.sh"
