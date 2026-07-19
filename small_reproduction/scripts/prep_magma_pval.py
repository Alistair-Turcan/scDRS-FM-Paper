#!/usr/bin/env python3
"""
Convert LDSC-format *.sumstats.gz -> MAGMA per-trait inputs.

For each trait, MAGMA gene analysis needs:
  - <trait>.pval : two columns  SNP  P   (P from the two-sided Z -> p)
  - <trait>.N    : the (median) sample size N used with `--pval ... N=<N>`

Input sumstats columns (LDSC): SNP A1 A2 N CHISQ Z  (Z is the signed z-score).
P = 2 * sf(|Z|) via the standard normal survival function.

Usage:
  python prep_magma_pval.py <sumstats_dir> <out_pval_dir>
Env override: none. Deterministic.
"""
import sys
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm

def main():
    sumstats_dir = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(sumstats_dir.glob("*.sumstats.gz"))
    print(f"found {len(files)} sumstats files")
    for f in files:
        trait = f.name.replace(".sumstats.gz", "")
        df = pd.read_csv(f, sep="\t")
        if "Z" not in df.columns or "SNP" not in df.columns:
            print(f"  [skip] {trait}: missing SNP/Z columns ({list(df.columns)})")
            continue
        z = df["Z"].astype(float).values
        p = 2.0 * norm.sf(np.abs(z))
        out = pd.DataFrame({"SNP": df["SNP"].values, "P": p})
        out.to_csv(out_dir / f"{trait}.pval", sep="\t", index=False)
        # N: use median (robust to per-SNP N variation); fall back to max
        if "N" in df.columns:
            N = int(round(float(np.nanmedian(df["N"].astype(float)))))
        else:
            N = 100000
        (out_dir / f"{trait}.N").write_text(str(N))
        print(f"  {trait}: {len(out):,} SNPs, N={N}")

if __name__ == "__main__":
    main()
