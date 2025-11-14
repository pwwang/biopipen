from contextlib import suppress
import hashlib
import shutil
import re
from pathlib import Path, PosixPath  # noqa: F401
from biopipen.utils.misc import run_command

fastqs: list[Path] = {{in.fastqs | each: as_path}}  # pyright: ignore  # noqa
outdir: Path = Path({{out.outdir | quote}})  # pyright: ignore
id: str = {{out.outdir | basename | quote}}  # pyright: ignore

cellranger = {{envs.cellranger | quote}}  # pyright: ignore
tmpdir = Path({{envs.tmpdir | quote}})  # pyright: ignore
ref: str = {{envs.ref | quote}}  # pyright: ignore
ncores = {{envs.ncores | int}}  # pyright: ignore
include_introns = {{envs.include_introns | repr}}  # pyright: ignore
create_bam = {{envs.create_bam | repr}}  # pyright: ignore
outdir_is_mounted: bool = {{envs.outdir_is_mounted | repr}}  # pyright: ignore
copy_outs_only: bool = {{envs.copy_outs_only | repr}}  # pyright: ignore

ref: Path = Path(ref).resolve()  # pyright: ignore
if not ref.exists():
    raise FileNotFoundError(f"Reference path does not exist: {ref}")
include_introns = str(include_introns).lower()
create_bam = str(create_bam).lower()

# create a temporary unique directory to store the soft-linked fastq files
uid = hashlib.md5(str(fastqs).encode()).hexdigest()[:8]
fastqdir = tmpdir / f"cellranger_count_{uid}"
fastqdir.mkdir(parents=True, exist_ok=True)
if len(fastqs) == 1 and fastqs[0].is_dir():
    fastqs = list(fastqs[0].glob("*.fastq.gz"))

# soft-link the fastq files to the temporary directory
for fastq in fastqs:
    fastq = Path(fastq)
    fqnames = re.split(r"(_S\d+_)", fastq.name)
    if len(fqnames) != 3:
        raise ValueError(
            fr"Expect one and only one '_S\d+_' in fastq file name: {fastq.name}"
        )

    linked = fastqdir / f"{id}{fqnames[1]}{fqnames[2]}"
    if linked.exists():
        linked.unlink()

    linked.symlink_to(fastq)

other_args = {{envs | dict_to_cli_args: dashify=True, exclude=['no_bam', 'create_bam', 'include_introns', 'cellranger', 'transcriptome', 'ref', 'tmpdir', 'id', 'ncores', 'outdir_is_mounted', 'copy_outs_only']}}  # pyright: ignore

command = [
    cellranger,
    "count",
    "--id",
    id,
    "--fastqs",
    fastqdir,
    "--transcriptome",
    str(ref),
    "--localcores",
    ncores,
    "--disable-ui",
    "--include-introns",
    include_introns,
    *other_args,
]

# check cellranger version
#   cellranger cellranger-7.2.0
version: str = run_command([cellranger, "--version"], stdout = "RETURN")  # type: ignore
version = version.replace("cellranger", "").replace("-", "").strip()  # type: ignore
print(f"# Detected cellranger version: {version}")
version: list[int] = list(map(int, version.split(".")))  # type: ignore
if version[0] >= 8:
    command += ["--create-bam", create_bam]
elif create_bam != "true":
    command += ["--no-bam"]

if outdir_is_mounted:
    print("# Using mounted outdir, redirecting cellranger output to a local tmpdir")
    local_outdir = tmpdir / f"{outdir.name}-{uid}" / id
    if local_outdir.parent.exists():
        shutil.rmtree(local_outdir.parent)
    local_outdir.parent.mkdir(parents=True, exist_ok=True)
    odir = local_outdir
else:
    odir = outdir

run_command(command, fg=True, cwd=str(odir.parent))

web_summary_html = odir / "outs" / "web_summary.html"
if not web_summary_html.exists():
    raise RuntimeError(
        f"web_summary.html does not exist in {odir}/outs. "
        "cellranger count failed."
    )

# Modify web_summary.html to move javascript to a separate file
# to void vscode live server breaking the page by injecting some code
print("# Modify web_summary.html to move javascript to a separate file")
try:
    web_summary_js = odir / "outs" / "web_summary.js"
    web_summary_content = web_summary_html.read_text()
    regex = re.compile(r"<script>(.+)</script>", re.DOTALL)
    web_summary_html.write_text(regex.sub(
        '<script src="web_summary.js"></script>',
        web_summary_content,
    ))
    web_summary_js.write_text(regex.search(web_summary_content).group(1))  # type: ignore
except Exception as e:
    print(f"Error modifying web_summary.html: {e}")
    raise e

# If using local tmpdir for output, move results to the final outdir
if outdir_is_mounted:
    print("# Copy results back to outdir")
    if outdir.exists():
        shutil.rmtree(outdir)

    if copy_outs_only:
        outdir.mkdir(parents=True, exist_ok=True)
        with suppress(Exception):
            # Some files may be failed to copy due to permission issues
            # But the contents are actually copied
            shutil.copytree(odir / "outs", outdir / "outs")
    else:
        with suppress(Exception):
            shutil.copytree(local_outdir, outdir)  # type: ignore

    # Make sure essential files exist
    web_summary_html = outdir / "outs" / "web_summary.html"
    web_summary_js = outdir / "outs" / "web_summary.js"
    for f in [web_summary_html, web_summary_js]:
        if not f.exists():
            raise RuntimeError(
                f"{f} does not exist in {outdir}/outs. "
                "Copying results back from tmpdir failed."
            )