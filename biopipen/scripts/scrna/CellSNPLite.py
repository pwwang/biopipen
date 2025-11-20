from __future__ import annotations

from contextlib import suppress
from pathlib import Path
from biopipen.core.filters import dict_to_cli_args
from biopipen.utils.misc import run_command

bamfile = Path({{in.bamfile | quote}})  # noqa: E999 # pyright: ignore
barcodefile = {{in.barcodefile | quote}}  # pyright: ignore
outdir = {{out.outdir | quote}}  # pyright: ignore
envs: dict = {{envs | repr}}  # pyright: ignore
cellsnp_lite = envs.pop("cellsnp_lite")
ncores = envs.pop("ncores")

with suppress(RuntimeError):
    run_command([cellsnp_lite, "--version"], fg=True)
    print("")

envs["nproc"] = ncores
envs["samFile"] = bamfile
envs["barcodeFile"] = barcodefile
envs["outDir"] = outdir

cmd = [cellsnp_lite, *dict_to_cli_args(envs)]
run_command(cmd, fg=True)
