#!/usr/bin/env python3
"""
Thin compatibility shim.

The original scDRS-FM SLURM scripts invoke `run_scdrs+.py`. The public package
ships the identical entry point as `run_scdrs_fm.py` (the '+' is not filesystem
friendly). This shim simply forwards all arguments to the UNMODIFIED package
entry point at scDRS-FM-main/run_scdrs_fm.py, so you can call either name.

It performs NO logic of its own and does not modify the package.

Resolution order for the package entry point:
  1. $SCDRS_FM_PY               (explicit path to run_scdrs_fm.py)
  2. $SCDRS_FM_HOME/run_scdrs_fm.py
  3. <repo>/scDRS-FM-main/run_scdrs_fm.py
"""
import os
import runpy
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
repo = here.parent

candidates = []
if os.environ.get("SCDRS_FM_PY"):
    candidates.append(Path(os.environ["SCDRS_FM_PY"]))
if os.environ.get("SCDRS_FM_HOME"):
    candidates.append(Path(os.environ["SCDRS_FM_HOME"]) / "run_scdrs_fm.py")
candidates.append(repo / "scDRS-FM-main" / "run_scdrs_fm.py")

target = next((c for c in candidates if c.is_file()), None)
if target is None:
    sys.stderr.write(
        "ERROR: could not locate the scDRS-FM package entry point run_scdrs_fm.py.\n"
        "Set SCDRS_FM_HOME to your scDRS-FM-main checkout, or place it under "
        f"{repo/'scDRS-FM-main'}.\nTried:\n  " + "\n  ".join(str(c) for c in candidates) + "\n"
    )
    sys.exit(2)

# Make the package importable and forward argv unchanged.
sys.path.insert(0, str(target.parent))
sys.argv = [str(target)] + sys.argv[1:]
runpy.run_path(str(target), run_name="__main__")
