from __future__ import annotations

from pathlib import Path
from contextlib import suppress
from biopipen.core.filters import dict_to_cli_args
from biopipen.utils.misc import run_command

cellsnpout = {{in.cellsnpout | quote}}  # noqa: E999 # pyright: ignore
outdir = {{out.outdir | quote}}  # pyright: ignore
envs: dict = {{envs | repr}}  # pyright: ignore
mquad = envs.pop("mquad")
ncores = envs.pop("ncores")
seed = envs.pop("seed", 8525)

with suppress(RuntimeError):
    run_command([mquad], fg=True)
    print("")

envs["cellData"] = cellsnpout
envs["outDir"] = outdir
envs["randSeed"] = seed
envs["nproc"] = ncores

cmd = [mquad, *dict_to_cli_args(envs, sep="=")]
run_command(cmd, fg=True)
