#!/usr/bin/env python3
"""Assemble per-trait MAGMA .genes.out files into the wide ZSTAT matrix that
scDRS-FM notebook 08 expects: MAGMA_v108_GENE_10_ZSTAT_for_scDRS.txt

Format (validated against build_gene_level_correlation in nb08):
  - TSV (sep='\t'), read with pd.read_csv(sep='\t')
  - one column per trait (column name == trait name)
  - rows = genes; values = MAGMA ZSTAT
  - gene index = HGNC symbol (scDRS convention), mapped from Entrez via NCBI37.3.gene.loc
The first column header is 'GENE' (index); nb08 selects only trait columns so the
index name is irrelevant to it, but symbols keep the file consistent with scDRS.
"""
import os
import sys
import pandas as pd
from pathlib import Path

OUT_DIR = Path(os.environ.get("MAGMA_OUT", os.path.join(os.environ.get("MAGMA_REF","./magma_ref"), "out")))
GENE_LOC = Path(os.environ.get("MAGMA_GENELOC", os.path.join(os.environ.get("MAGMA_REF","./magma_ref"), "gene_loc/NCBI37.3.gene.loc")))
DEST_LOCAL = Path(os.environ.get("MAGMA_ZSTAT_OUT", os.path.join(os.environ.get("MAGMA_REF","./magma_ref"), "MAGMA_v108_GENE_10_ZSTAT_for_scDRS.txt")))
DEST_DURABLE = DEST_LOCAL  # single destination in the repo layout

def main():
    # Entrez -> symbol map (col0=Entrez, col5=symbol; no header)
    loc = pd.read_csv(GENE_LOC, sep="\t", header=None,
                      usecols=[0, 5], names=["entrez", "symbol"], dtype={0: str})
    ent2sym = dict(zip(loc["entrez"].astype(str), loc["symbol"]))
    print(f"gene.loc: {len(ent2sym)} Entrez->symbol mappings")

    files = sorted(OUT_DIR.glob("*.genes.out"))
    print(f"found {len(files)} .genes.out files")
    series = {}
    for f in files:
        trait = f.name.replace(".genes.out", "")
        # MAGMA .genes.out is whitespace-delimited
        df = pd.read_csv(f, sep=r"\s+")
        df["GENE"] = df["GENE"].astype(str)
        df["SYMBOL"] = df["GENE"].map(ent2sym)
        df = df.dropna(subset=["SYMBOL"])
        # collapse duplicate symbols by max |ZSTAT| (keep most significant)
        df = df.reindex(df["ZSTAT"].abs().sort_values(ascending=False).index)
        df = df.drop_duplicates(subset="SYMBOL", keep="first")
        series[trait] = df.set_index("SYMBOL")["ZSTAT"]
    mat = pd.DataFrame(series).sort_index()
    mat.index.name = "GENE"
    print(f"matrix: {mat.shape[0]} genes x {mat.shape[1]} traits")
    print(f"non-null cells: {int(mat.notna().sum().sum())} / {mat.size}")
    # write with GENE as first column
    mat.to_csv(DEST_LOCAL, sep="\t")
    print(f"wrote {DEST_LOCAL}")
    # copy to durable (S3-FUSE: copy works, not rename)
    import shutil
    if DEST_DURABLE != DEST_LOCAL:
        shutil.copyfile(DEST_LOCAL, DEST_DURABLE)
        print(f"copied to {DEST_DURABLE}")
    # quick sanity: show a few known immune traits present
    for t in ["PASS_Celiac", "PASS_Rheumatoid_Arthritis", "PASS_Alzheimers_Jansen2019"]:
        if t in mat.columns:
            print(f"  {t}: {mat[t].notna().sum()} genes with ZSTAT")

if __name__ == "__main__":
    main()
