import hashlib
import shutil
import shlex
import re
from contextlib import suppress
from pathlib import Path, PosixPath  # noqa: F401
from panpath import LocalPath  # noqa: F401
from biopipen.utils.misc import run_command

csv: Path = Path({{in.csv | quote}})  # pyright: ignore
outdir: Path = Path({{out.outdir | quote}})  # pyright: ignore
id: str = {{out.outdir | basename | quote}}  # pyright: ignore

cellranger: str = {{envs.cellranger | quote}}  # pyright: ignore
tmpdir = Path({{envs.tmpdir | quote}})  # pyright: ignore
ncores: int = {{envs.ncores | int}}  # pyright: ignore
outdir_is_mounted: bool = {{envs.outdir_is_mounted | repr}}  # pyright: ignore
copy_outs_only: bool = {{envs.copy_outs_only | repr}}  # pyright: ignore

other_args = {{envs | dict_to_cli_args: dashify=True, exclude=['cellranger', 'tmpdir', 'ncores', 'outdir_is_mounted', 'copy_outs_only']}}  # pyright: ignore

uid = hashlib.md5(str(csv).encode()).hexdigest()[:8]

command = [
    *shlex.split(cellranger),
    "multi",
    "--id",
    id,
    "--csv",
    str(csv.resolve()),
    "--localcores",
    ncores,
    "--disable-ui",
    "--nopreflight",
    *other_args,
]

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

per_sample_outs = odir / "outs" / "per_sample_outs"
if not per_sample_outs.exists() or not any(per_sample_outs.iterdir()):
    raise RuntimeError(
        f"outs/per_sample_outs does not exist or is empty in {odir}. "
        "cellranger multi failed."
    )

# Modify web_summary.html files (per-sample and run-level) to move javascript
# to a separate file to void vscode live server breaking the page by injecting some code
print("# Modify web_summary.html files to move javascript to separate files")
html_files = list(per_sample_outs.glob("*/web_summary.html"))
qc_report = odir / "outs" / "qc_report.html"
if qc_report.exists():
    html_files.append(qc_report)

for html_file in html_files:
    try:
        js_file = html_file.with_name(html_file.stem + ".js")
        content = html_file.read_text()
        regex = re.compile(r"<script>(.+)</script>", re.DOTALL)
        m = regex.search(content)
        if m:
            html_file.write_text(regex.sub(
                f'<script src="{js_file.name}"></script>',
                content,
            ))
            js_file.write_text(m.group(1))
    except Exception as e:
        print(f"Warning: Error modifying {html_file}: {e}")

# If using local tmpdir for output, move results to the final outdir
if outdir_is_mounted:
    print("# Copy results back to outdir")
    if outdir.exists():
        shutil.rmtree(outdir)

    if copy_outs_only:
        outdir.mkdir(parents=True, exist_ok=True)
        with suppress(Exception):
            shutil.copytree(odir / "outs", outdir / "outs")
    else:
        with suppress(Exception):
            shutil.copytree(local_outdir, outdir)  # type: ignore

    # Verify per_sample_outs was copied
    per_sample_outs = outdir / "outs" / "per_sample_outs"
    if not per_sample_outs.exists() or not any(per_sample_outs.iterdir()):
        raise RuntimeError(
            f"outs/per_sample_outs does not exist or is empty in {outdir}. "
            "Copying results back from tmpdir failed."
        )
