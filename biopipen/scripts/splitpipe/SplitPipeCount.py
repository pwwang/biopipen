import shlex
import subprocess
import sys
from pathlib import Path, PosixPath  # noqa: F401
from panpath import LocalPath  # noqa: F401

fq1: list[Path] = {{in.fq1 | each: as_path}}  # pyright: ignore  # noqa
fq2: list[Path] = {{in.fq2 | default: [] | each: as_path}}  # pyright: ignore  # noqa
outdir: Path = Path({{out.outdir | quote}})  # pyright: ignore

splitpipe = {{envs.splitpipe | quote}}  # pyright: ignore
nthreads = {{ (envs.ncores if envs.ncores is not none else envs.nthreads) | int }}  # pyright: ignore
chemistry = {{envs.chemistry | repr}}  # pyright: ignore
kit = {{envs.kit | repr}}  # pyright: ignore
genome_dir = {{in.genome_dir | default: None | quote}} or {{envs.genome_dir | quote}}  # pyright: ignore
samples = {{envs.samples | repr}}  # pyright: ignore
samp_sltab = {{envs.samp_sltab | repr}}  # pyright: ignore
samp_slt = {{envs.samp_slt | repr}}  # pyright: ignore
samp_list = {{envs.samp_list | repr}}  # pyright: ignore
parfile = {{envs.parfile | repr}}  # pyright: ignore
targeted_list = {{envs.targeted_list | repr}}  # pyright: ignore
crispr: bool = {{envs.crispr | repr}}  # pyright: ignore
crsp_guides = {{envs.crsp_guides | repr}}  # pyright: ignore
crsp_read_thresh = {{envs.crsp_read_thresh | repr}}  # pyright: ignore
crsp_tscp_thresh = {{envs.crsp_tscp_thresh | repr}}  # pyright: ignore
ffpe: bool = {{envs.ffpe | repr}}  # pyright: ignore
tcr_analysis: bool = {{envs.tcr_analysis | repr}}  # pyright: ignore
bcr_analysis: bool = {{envs.bcr_analysis | repr}}  # pyright: ignore
immune_genome = {{envs.immune_genome | repr}}  # pyright: ignore
parent_dir = {{envs.parent_dir | repr}}  # pyright: ignore
use_imgt_db: bool = {{envs.use_imgt_db | repr}}  # pyright: ignore
min_map_fraction = {{envs.min_map_fraction | repr}}  # pyright: ignore

other_args = {{envs | dict_to_cli_args: exclude=['splitpipe', 'nthreads', 'ncores', 'chemistry', 'kit', 'genome_dir', 'samples', 'samp_sltab', 'samp_slt', 'samp_list', 'parfile', 'targeted_list', 'crispr', 'crsp_guides', 'crsp_read_thresh', 'crsp_tscp_thresh', 'ffpe', 'tcr_analysis', 'bcr_analysis', 'immune_genome', 'parent_dir', 'use_imgt_db', 'min_map_fraction']}}  # pyright: ignore

if not chemistry:
    sys.exit("`envs.chemistry` is required by `split-pipe --mode all`.")
if not kit:
    sys.exit("`envs.kit` is required by `split-pipe --mode all`.")
if not genome_dir:
    sys.exit(
        "`envs.genome_dir` (or the `genome_dir` input) is required by "
        "`split-pipe --mode all`."
    )
if not fq1:
    sys.exit("No fastq (R1) files provided for the sublibrary.")
if not Path(genome_dir).exists():
    sys.exit(f"Genome directory does not exist: {genome_dir}")

for f in [*fq1, *fq2]:
    if not Path(f).exists():
        sys.exit(f"Input fastq file does not exist: {f}")

if crispr and not (crsp_guides and parent_dir):
    sys.exit(
        "`envs.crsp_guides` and `envs.parent_dir` are required when "
        "`envs.crispr` is True."
    )

if samp_sltab and samp_slt:
    sys.exit("`envs.samp_sltab` and `envs.samp_slt` cannot be used together.")

if min_map_fraction is not None and parfile:
    sys.exit(
        "`envs.min_map_fraction` and `envs.parfile` cannot be used together, "
        "since both define `post_min_map_frac`. Put `post_min_map_frac` in "
        "`envs.parfile` instead, or unset `envs.parfile`."
    )

# Parse documents the minimum mapping fraction as `post_min_map_frac` in the
# parameter file, which is the only documented way to set it
if min_map_fraction is not None:
    parfile = outdir.parent / f"{outdir.name}-post_min_map_frac.parfile"
    parfile.write_text(f"post_min_map_frac {min_map_fraction}\n")
    print(f"# Wrote parameter file: {parfile}")

command = [
    splitpipe,
    "--mode",
    "all",
    "--chemistry",
    chemistry,
    "--kit",
    kit,
    "--nthreads",
    str(nthreads),
    "--fq1",
    *[str(f) for f in fq1],
    "--output_dir",
    str(outdir),
    "--genome_dir",
    str(genome_dir),
]

if fq2:
    command += ["--fq2", *[str(f) for f in fq2]]

if samples:
    for sample in samples:
        fields = str(sample).split()
        if len(fields) != 2:
            sys.exit(
                "Each entry of `envs.samples` should be `NAME WELLS`, "
                f"e.g. `xcond_1 A1-A6`. Got: {sample}"
            )
        command += ["--sample", *fields]

if samp_sltab:
    command += ["--samp_sltab", str(samp_sltab)]

if samp_slt:
    command += ["--samp_slt", str(samp_slt)]

if samp_list:
    command += ["--samp_list", str(samp_list)]

if parfile:
    command += ["--parfile", str(parfile)]

if targeted_list:
    command += ["--targeted_list", str(targeted_list)]

if crispr:
    command += ["--crispr", "--crsp_guides", str(crsp_guides)]

if crsp_read_thresh is not None:
    command += ["--crsp_read_thresh", str(crsp_read_thresh)]

if crsp_tscp_thresh is not None:
    command += ["--crsp_tscp_thresh", str(crsp_tscp_thresh)]

if ffpe:
    command += ["--ffpe"]

if tcr_analysis:
    command += ["--tcr_analysis"]

if bcr_analysis:
    command += ["--bcr_analysis"]

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
        f"# split-pipe --mode all exited with code {proc.returncode}. "
        f"stderr tail:\n{tail}",
        file=sys.stderr,
    )
    sys.exit(proc.returncode)

print("\n".join(stdout_file.read_text().strip().splitlines()[-30:]))
