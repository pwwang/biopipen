import shlex
import subprocess
import sys
from pathlib import Path, PosixPath  # noqa: F401
from panpath import LocalPath  # noqa: F401

fasta: list[Path] = {{in.fasta | each: as_path}}  # pyright: ignore  # noqa
genes: list[Path] = {{in.genes | each: as_path}}  # pyright: ignore  # noqa
outdir: Path = Path({{out.outdir | quote}})  # pyright: ignore

splitpipe: str = {{envs.splitpipe | quote}}  # pyright: ignore
nthreads: int = {{ (envs.ncores if envs.ncores is not none else envs.nthreads) | int }}  # pyright: ignore
genome_name: str = {{envs.genome_name | repr}}  # pyright: ignore

other_args = {{envs | dict_to_cli_args: exclude=['splitpipe', 'nthreads', 'ncores', 'genome_name']}}  # pyright: ignore

if not genome_name:
    sys.exit("`envs.genome_name` is required to build a reference genome.")

genome_names: list[str] = (
    [genome_name] if isinstance(genome_name, str) else list(genome_name)
)
if len(genome_names) != len(fasta) or len(fasta) != len(genes):
    sys.exit(
        "The number of genome names, fasta files and gene annotation files "
        "must match: "
        f"{len(genome_names)} genome name(s), {len(fasta)} fasta file(s), "
        f"{len(genes)} gene annotation file(s)."
    )

for f in [*fasta, *genes]:
    if not Path(f).exists():
        sys.exit(f"Input file does not exist: {f}")

command = [
    splitpipe,
    "--mode",
    "mkref",
    "--genome_name",
    *genome_names,
    "--nthreads",
    str(nthreads),
    "--fasta",
    *[str(f) for f in fasta],
    "--genes",
    *[str(f) for f in genes],
    "--output_dir",
    str(outdir),
    *other_args,
]

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
        f"# split-pipe --mode mkref exited with code {proc.returncode}. "
        f"stderr tail:\n{tail}",
        file=sys.stderr,
    )
    sys.exit(proc.returncode)

print("\n".join(stdout_file.read_text().strip().splitlines()[-30:]))
