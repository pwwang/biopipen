from contextlib import suppress
import hashlib
import re
from pathlib import Path, PosixPath  # noqa: F401
from panpath import LocalPath  # noqa: F401
from biopipen.utils.misc import run_command

fastqs: list[Path] = {{in.fastqs | each: as_path}}  # pyright: ignore  # noqa
outdir: Path = Path({{out.outdir | quote}})  # pyright: ignore
id: str = {{out.outdir | basename | quote}}  # pyright: ignore

ncores = {{envs.ncores | int}}  # pyright: ignore
pipseeker = {{envs.pipseeker | quote}}  # pyright: ignore
ref: str = {{envs.ref | quote}}  # pyright: ignore
verbosity = {{envs.verbosity | int}}  # pyright: ignore
remove_bam = {{envs.remove_bam | repr}}  # pyright: ignore
skip_version_check = {{envs.skip_version_check | repr}}  # pyright: ignore
chemistry = {{envs.chemistry | quote}}  # pyright: ignore
tmpdir = Path({{envs.tmpdir | quote}})  # pyright: ignore

ref: Path = Path(ref).resolve()  # pyright: ignore
if not ref.exists():
    raise FileNotFoundError(f"Reference path does not exist: {ref}")

# create a temporary unique directory to store the soft-linked fastq files
uid = hashlib.md5(str(fastqs).encode()).hexdigest()[:8]
fastqdir = tmpdir / f"pipseeker_full_{uid}"
fastqdir.mkdir(parents=True, exist_ok=True)
if len(fastqs) == 1 and fastqs[0].is_dir():
    fastqs = list(fastqs[0].glob("*.fastq.gz"))

# soft-link the fastq files to the temporary directory
for fastq in fastqs:
    fastq = Path(fastq)
    fqnames = re.split(r"(_R\d+)", fastq.name)
    if len(fqnames) != 3:
        raise ValueError(
            fr"Expect one and only one '_R\d+' in fastq file name: {fastq.name}"
        )

    linked = fastqdir / f"{id}{fqnames[1]}{fqnames[2]}"
    if linked.exists():
        linked.unlink()

    linked.symlink_to(fastq)

other_args = {{envs | dict_to_cli_args: dashify=True, exclude=['ncores', 'pipseeker', 'ref', 'verbosity', 'remove_bam', 'skip_version_check', 'chemistry', 'tmpdir']}}  # pyright: ignore

command = [
    pipseeker,
    "full",
    "--fastq",
    f"{fastqdir}/{id}",
    "--chemistry",
    chemistry,
    "--threads",
    ncores,
    "--star-index-path",
    str(ref),
    "--verbosity",
    verbosity,
    "--output-path",
    outdir,
    *other_args,
]

run_command(command, fg=True)
