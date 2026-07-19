#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# download_figshare_data.sh
#
# Downloads the scDRS-FM data release from figshare (article 33000602):
#   - 8 dataset bundles, each an .h5ad expression matrix + matching .cov file
#       TMS_FACS, TS_FACS, TMS_Droplet, SEA_AD, Soskic, Nathan, Cano_Gamez, Braun
#   - GS_files_75.zip                    -> the 75 MAGMA gene sets (= gs_split/)
#   - Microglia_phenotype_genesets.zip   -> 5 microglia signatures
#   - T_cell_phenotype_genesets.zip      -> T-cell phenotype signatures
#   - genet_cor.csv                      -> LDSC genetic-correlation (rg) matrix
#   - 75 *.sumstats.gz                   -> LDSC-format GWAS summary statistics
#
# The full release is ~10.2 GB compressed (~26 GB uncompressed h5ads).
#
# USAGE:
#   bash download/download_figshare_data.sh [DEST_DIR]
# DEST_DIR defaults to ./data
#
# NOTE: figshare download URLs are resolved from the article's public API.
# If the article layout changes, edit FIGSHARE_ARTICLE below or download the
# files manually from https://figshare.com/articles/dataset/33000602 and place
# them under DEST_DIR following download/manifest.tsv.
# ---------------------------------------------------------------------------
set -euo pipefail

DEST="${1:-./data}"
FIGSHARE_ARTICLE=33000602
API="https://api.figshare.com/v2/articles/${FIGSHARE_ARTICLE}"

mkdir -p "${DEST}/raw" "${DEST}/extracted" "${DEST}/sumstats"

echo "[figshare] Querying article ${FIGSHARE_ARTICLE} file list..."
# Requires: curl, jq (jq optional — falls back to python if absent)
if command -v jq >/dev/null 2>&1; then
  curl -sSL "${API}" | jq -r '.files[] | "\(.download_url)\t\(.name)\t\(.size)"' > "${DEST}/raw/_filelist.tsv"
else
  curl -sSL "${API}" | python3 -c '
import sys, json
d = json.load(sys.stdin)
for f in d["files"]:
    print(f"{f[\"download_url\"]}\t{f[\"name\"]}\t{f[\"size\"]}")
' > "${DEST}/raw/_filelist.tsv"
fi

echo "[figshare] Files advertised by the article:"
cut -f2,3 "${DEST}/raw/_filelist.tsv" | sed 's/^/    /'

echo "[figshare] Downloading (resumable)..."
while IFS=$'\t' read -r url name size; do
  out="${DEST}/raw/${name}"
  if [[ -f "${out}" ]]; then
    echo "  [skip] ${name} already present"
    continue
  fi
  echo "  [get ] ${name} (${size} bytes)"
  curl -sSL -C - -o "${out}" "${url}"
done < "${DEST}/raw/_filelist.tsv"

echo "[figshare] Unzipping dataset bundles + gene-set archives..."
shopt -s nullglob
for z in "${DEST}"/raw/*.zip; do
  echo "  [unzip] $(basename "$z")"
  unzip -n -q "$z" -d "${DEST}/extracted"
done

# Gene sets: GS_files_75 -> gs_split ; phenotype archives -> gene_sets/
mkdir -p "${DEST}/gene_sets"
if [[ -d "${DEST}/extracted/gs_split" ]]; then
  ln -sfn "${DEST}/extracted/gs_split" "${DEST}/gene_sets/gs_split"
fi

# Move any *.sumstats.gz into a single dir
find "${DEST}/raw" "${DEST}/extracted" -maxdepth 2 -name "*.sumstats.gz" -exec cp -n {} "${DEST}/sumstats/" \; 2>/dev/null || true

echo "[figshare] Done. Verify contents against download/manifest.tsv:"
echo "    Expected: 8 datasets (h5ad + cov), GS_files_75, 2 phenotype archives, genet_cor.csv, 75 sumstats.gz"
echo "[figshare] Extracted tree:"
ls -la "${DEST}/extracted" | sed 's/^/    /'
