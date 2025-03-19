import uuid
import re
import os.path
from pathlib import Path, PosixPath  # noqa: F401
from biopipen.utils.misc import run_command

fastqs: list[Path] = {{in.fastqs | each: as_path}}  # pyright: ignore  # noqa
outdir: str = {{out.outdir | quote}}  # pyright: ignore
id = {{out.outdir | basename | quote}}  # pyright: ignore

cellranger = {{envs.cellranger | quote}}  # pyright: ignore
tmpdir = Path({{envs.tmpdir | quote}})  # pyright: ignore
ref: str = {{envs.ref | quote}}  # pyright: ignore
ncores = {{envs.ncores | int}}  # pyright: ignore
include_introns = {{envs.include_introns | repr}}  # pyright: ignore
create_bam = {{envs.create_bam | repr}}  # pyright: ignore

include_introns = str(include_introns).lower()
create_bam = str(create_bam).lower()

# create a temporary unique directory to store the soft-linked fastq files
fastqdir = tmpdir / f"cellranger_count_{uuid.uuid4()}"
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

other_args = {{envs | dict_to_cli_args: dashify=True, exclude=['no_bam', 'create_bam', 'include_introns', 'cellranger', 'transcriptome', 'ref', 'tmpdir', 'id', 'ncores']}}  # pyright: ignore

command = [
    cellranger,
    "count",
    "--id",
    id,
    "--fastqs",
    fastqdir,
    "--transcriptome",
    Path(ref).resolve(),
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
version: list[int] = list(map(int, version.split(".")))  # type: ignore
if version[0] >= 8:
    command += ["--create-bam", create_bam]
elif create_bam != "true":
    command += ["--no-bam"]

run_command(command, fg=True, cwd=str(Path(outdir).parent))

web_summary_html = Path(outdir) / "outs" / "web_summary.html"
if not web_summary_html.exists():
    raise RuntimeError(
        f"web_summary.html does not exist in {outdir}/outs. "
        "cellranger count failed."
    )

# Modify web_summary.html to move javascript to a separate file
# to void vscode live server breaking the page by injecting some code
print("# Modify web_summary.html to move javascript to a separate file")
try:
    web_summary_js = Path(outdir) / "outs" / "web_summary.js"
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
