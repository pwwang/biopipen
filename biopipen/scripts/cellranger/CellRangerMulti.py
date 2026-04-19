import hashlib
import shutil
import shlex
import re
from contextlib import suppress
from diot import Diot  # type: ignore
from pathlib import Path, PosixPath  # noqa: F401
from panpath import LocalPath  # noqa: F401
from biopipen.utils.misc import run_command

fastqs: list[Path] = {{in.fastqs | each: as_path}}  # pyright: ignore  # noqa
outdir: Path = Path({{out.outdir | quote}})  # pyright: ignore
id: str = {{out.outdir | basename | quote}}  # pyright: ignore

cellranger: str = {{envs.cellranger | quote}}  # pyright: ignore
tmpdir = Path({{envs.tmpdir | quote}})  # pyright: ignore
ncores: int = {{envs.ncores | int}}  # pyright: ignore
outdir_is_mounted: bool = {{envs.outdir_is_mounted | repr}}  # pyright: ignore
copy_outs_only: bool = {{envs.copy_outs_only | repr}}  # pyright: ignore
gex = {{envs.gex | repr}}  # pyright: ignore
vdj = {{envs.vdj | repr}}  # pyright: ignore
feature = {{envs.feature | repr}}  # pyright: ignore
libraries = {{envs.libraries | repr}}  # pyright: ignore

uid = hashlib.md5(str(fastqs).encode()).hexdigest()[:8]

# ── Symlink all input FASTQs into a temporary directory ──────────────────────
fastqdir = tmpdir / f"cellranger_multi_{uid}"
fastqdir.mkdir(parents=True, exist_ok=True)

if len(fastqs) == 1 and fastqs[0].is_dir():
    fastqs = list(fastqs[0].glob("*.fastq.gz"))

for fastq in fastqs:
    linked = fastqdir / fastq.name
    if not linked.exists():
        linked.symlink_to(fastq.resolve())

# ── Build multi config CSV ────────────────────────────────────────────────────
csv_path = tmpdir / f"cellranger_multi_{uid}.csv"

def _write_section(f, section_name, section_dict):
    """Write a config section to the CSV file."""
    f.write(f"[{section_name}]\n")
    for key, val in section_dict.items():
        csv_key = key.replace("_", "-")
        if val is None or val is False:
            continue
        if val is True:
            val = "true"
        f.write(f"{csv_key},{val}\n")
    f.write("\n")

with open(csv_path, "w") as f:
    # [gene-expression] section
    if gex:
        gex_dict = dict(gex)
        ref = gex_dict.get("reference") or gex_dict.get("ref") or ""
        probe_set = gex_dict.get("probe_set") or gex_dict.get("probe-set") or ""
        if not ref and not probe_set:
            raise ValueError(
                "Reference genome not configured for [gene-expression]. "
                "Set 'envs.gex.reference' in your CellRangerMulti configuration, "
                "or set 'ref.ref_cellranger_gex' in your biopipen configuration "
                "(~/.biopipen.toml or ./.biopipen.toml)."
            )
        _write_section(f, "gene-expression", gex_dict)

    # [vdj] section
    if vdj:
        _write_section(f, "vdj", dict(vdj))

    # [feature] section
    if feature:
        _write_section(f, "feature", dict(feature))

    # [libraries] section
    has_lanes = any("lanes" in lib for lib in libraries)
    f.write("[libraries]\n")
    if has_lanes:
        f.write("fastq_id,fastqs,lanes,feature_types\n")
    else:
        f.write("fastq_id,fastqs,feature_types\n")
    for lib in libraries:
        lib_fastqs = lib.get("fastqs") or str(fastqdir)
        fastq_id = lib["fastq_id"]
        feature_types = lib["feature_types"]
        if has_lanes:
            lanes = lib.get("lanes", "")
            f.write(f"{fastq_id},{lib_fastqs},{lanes},{feature_types}\n")
        else:
            f.write(f"{fastq_id},{lib_fastqs},{feature_types}\n")

print(f"# Generated multi config CSV at: {csv_path}")
print(open(csv_path).read())

command = [
    *shlex.split(cellranger),
    "multi",
    "--id",
    id,
    "--csv",
    str(csv_path),
    "--localcores",
    ncores,
    "--disable-ui",
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
