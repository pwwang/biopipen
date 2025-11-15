import hashlib
import shutil
import re
from contextlib import suppress
from pathlib import Path, PosixPath  # noqa: F401
from biopipen.utils.misc import run_command

fastqs: list[Path] = {{in.fastqs | each: as_path}}  # pyright: ignore  # noqa
outdir: Path = Path({{out.outdir | quote}})  # pyright: ignore
id: str = {{out.outdir | basename | quote}}  # pyright: ignore

cellranger: str = {{envs.cellranger | quote}}  # pyright: ignore
tmpdir = Path({{envs.tmpdir | quote}})  # pyright: ignore
ref: str = {{envs.ref | quote}}  # pyright: ignore
ncores: int = {{envs.ncores | int}}  # pyright: ignore
outdir_is_mounted: bool = {{envs.outdir_is_mounted | repr}}  # pyright: ignore
copy_outs_only: bool = {{envs.copy_outs_only | repr}}  # pyright: ignore

# create a temporary unique directory to store the soft-linked fastq files
uid = hashlib.md5(str(fastqs).encode()).hexdigest()[:8]
fastqdir = tmpdir / f"cellranger_count_{uid}"
fastqdir.mkdir(parents=True, exist_ok=True)
if len(fastqs) == 1 and fastqs[0].is_dir():
    fastqs = list(fastqs[0].glob("*.fastq.gz"))

# soft-link the fastq files to the temporary directory
for fastq in fastqs:
    fastq = Path(fastq)
    (fastqdir / fastq.name).symlink_to(fastq)

other_args = {{envs | dict_to_cli_args: dashify=True, exclude=['cellranger', 'reference', 'ref', 'tmpdir', 'id', 'ncores', 'outdir_is_mounted', 'copy_outs_only']}}  # pyright: ignore

command = [
    cellranger,
    "vdj",
    "--id",
    id,
    "--fastqs",
    fastqdir,
    "--reference",
    Path(ref).resolve(),
    "--localcores",
    ncores,
    "--disable-ui",
    *other_args,
]

version: str = run_command([cellranger, "--version"], stdout = "RETURN")  # type: ignore
version = version.replace("cellranger", "").replace("-", "").strip()  # type: ignore
print(f"# Detected cellranger version: {version}")

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
        "cellranger vdj failed."
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
    filtered_annotations_csv = outdir / "outs" / "filtered_contig_annotations.csv"
    for f in [web_summary_html, web_summary_js, filtered_annotations_csv]:
        if not f.exists():
            raise RuntimeError(
                f"{f} does not exist in {outdir}/outs. "
                "Copying results back from tmpdir failed."
            )
