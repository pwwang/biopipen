from __future__ import annotations

import shutil
import sys
from pathlib import Path
from contextlib import suppress
from biopipen.core.filters import dict_to_cli_args
from biopipen.utils.misc import run_command
from biopipen import __file__ as biopipen_file

cellsnpout = {{in.cellsnpout | quote}}  # noqa: E999 # pyright: ignore
outdir = {{out.outdir | quote}}  # pyright: ignore
envs: dict = {{envs | repr}}  # pyright: ignore
mquad: str = envs.pop("mquad") or "mquad"
ncores = envs.pop("ncores")
seed = envs.pop("seed", 8525)
mquad_cli = Path(biopipen_file).parent / "scripts" / "scrna" / "mquad_cli.py"

with suppress(RuntimeError):
    run_command([mquad], fg=True)
    print("")

# Patch mquad to fix the following issue:
# Traceback (most recent call last):
#   File ".../python3.12/multiprocessing/pool.py", line 125, in worker
#     result = (True, func(*args, **kwds))
#                     ^^^^^^^^^^^^^^^^^^^
#   File ".../python3.12/multiprocessing/pool.py", line 51, in starmapstar
#     return list(itertools.starmap(args[0], args[1]))
#            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File ".../python3.12/site-packages/mquad/mquad_batch_mixbin.py", line 307, in fit_batch
#     basic_stats = basicStats(valid_ad, valid_dp, valid_row_sizes)
#                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File ".../python3.12/site-packages/mquad/mquad_batch_mixbin.py", line 321, in basicStats
#     _d = valid_dp[left:right]
#          ~~~~~~~~^^^^^^^^^^^^
# IndexError: too many indices for array: array is 0-dimensional, but 1 were indexed

# get python path from mquad executable
mquad = mquad if "/" in mquad else shutil.which(mquad)  # type: ignore
if mquad is None:
    raise FileNotFoundError("mquad executable not found in PATH")

py_path = sys.executable
with open(mquad, "r") as f:
    for line in f:
        if line.startswith("#!"):
            py_path = line[2:].strip()
            break

envs["cellData"] = cellsnpout
envs["outDir"] = outdir
envs["randSeed"] = seed
envs["nproc"] = ncores

cmd = [py_path, mquad_cli, *dict_to_cli_args(envs, sep="=")]
run_command(cmd, fg=True, bufsize=1)
