import shlex
import subprocess
import sys
from pathlib import Path, PosixPath  # noqa: F401
from panpath import LocalPath  # noqa: F401

sublibs: list[Path] = {{in.sublibs | each: as_path}}  # pyright: ignore  # noqa
outdir: Path = Path({{out.outdir | quote}})  # pyright: ignore

splitpipe = {{envs.splitpipe | quote}}  # pyright: ignore
nthreads = {{ (envs.ncores if envs.ncores is not none else envs.nthreads) | int }}  # pyright: ignore
sublib_list = {{envs.sublib_list | repr}}  # pyright: ignore
immune_genome = {{envs.immune_genome | repr}}  # pyright: ignore
parent_dir = {{envs.parent_dir | repr}}  # pyright: ignore
use_imgt_db: bool = {{envs.use_imgt_db | repr}}  # pyright: ignore

other_args = {{envs | dict_to_cli_args: exclude=['splitpipe', 'nthreads', 'ncores', 'sublib_list', 'immune_genome', 'parent_dir', 'use_imgt_db']}}  # pyright: ignore

if not sublib_list and not sublibs:
    sys.exit(
        "No sublibrary to combine. Provide `in.sublibs` (or `envs.sublib_list`)."
    )

for sublib in sublibs:
    if not Path(sublib).exists():
        sys.exit(f"Sublibrary directory does not exist: {sublib}")

command = [
    splitpipe,
    "--mode",
    "comb",
    "--nthreads",
    str(nthreads),
    "--output_dir",
    str(outdir),
]

if sublib_list:
    command += ["--sublib_list", str(sublib_list)]
else:
    command += ["--sublibraries", *[str(s) for s in sublibs]]

if immune_genome:
    command += ["--immune_genome", immune_genome]

if parent_dir:
    command += ["--parent_dir", str(parent_dir)]

if use_imgt_db:
    command += ["--use_imgt_db"]

command += other_args

print("# Running: " + " ".join(shlex.quote(c) for c in command))
# split-pipe spawns a memory-profiler daemon (mem_loop.py) that inherits
# split-pipe's stdout/stderr and keeps them open. capture_output=True would
# block reading those pipes until EOF, which never arrives (the daemon
# outlives split-pipe). Redirect to files instead — do not revert this.
stdout_file = outdir / "split-pipe.stdout.txt"
stderr_file = outdir / "split-pipe.stderr.txt"
with stdout_file.open("w") as fout, stderr_file.open("w") as ferr:
    proc = subprocess.run(command, text=True, stdout=fout, stderr=ferr, check=False)
if proc.returncode != 0:
    tail = "\n".join(stderr_file.read_text().strip().splitlines()[-30:])
    print(
        f"# split-pipe --mode comb exited with code {proc.returncode}. "
        f"stderr tail:\n{tail}",
        file=sys.stderr,
    )
    sys.exit(proc.returncode)

print("\n".join(stdout_file.read_text().strip().splitlines()[-30:]))
