from pathlib import Path, PosixPath  # noqa: F401

from biopipen.utils.misc import logger
from biopipen.scripts.vcf.bcftools_utils import run_bcftools

infile: str | Path = {{in.infile | quote}}  # pyright: ignore # noqa: #999
outfile: str = {{out.outfile | quote}}  # pyright: ignore
outdir = Path(outfile).parent

envs: dict = {{envs | dict | repr}}  # pyright: ignore
bcftools = envs.pop("bcftools")
tabix = envs.pop("tabix")
keep = envs.pop("keep")
ncores = envs.pop("ncores")
includes = envs.pop("includes")
excludes = envs.pop("excludes")
gz = envs.pop("gz")
index = envs.pop("index")

# a.vcf.gz -> a
# a.vcf -> a
stem = Path(infile).stem
if stem.endswith(".vcf"):
    stem = stem[:-4]
# .vcf.gz
# .gz
ext = ".vcf.gz" if index or gz else '.vcf'


def normalize_expr(expr, flag, prev_n_filters=0):
    out = {}
    if not expr:
        return out
    if isinstance(expr, list):
        for ex in expr:
            out[f"FILTER_{flag.upper()}_{len(out) + 1 + prev_n_filters}"] = (ex, flag)
    elif isinstance(expr, dict):
        for name, ex in expr.items():
            out[name] = (ex, flag)
    else: # str
        out[f"FILTER_{flag.upper()}_{len(out) + 1 + prev_n_filters}"] = (expr, flag)
    return out


def handle_filter(vcf, fname, filt, flag, final):
    logger.info("- Handling filter %s: %s ...", fname, filt)

    arguments = envs.copy()
    arguments[flag] = filt
    arguments["_"] = vcf
    arguments["o"] = outfile if final else outdir / f"{stem}.{fname}{ext}"
    if keep:
        arguments["s"] = fname

    run_bcftools(arguments, bcftools=bcftools, index=index and final, tabix=tabix)

    if final:
        flagfile = outdir.joinpath(f"{stem}.{fname}{ext}")
        if flagfile.is_symlink():
            flagfile.unlink()
        outdir.joinpath(f"{stem}.{fname}{ext}").symlink_to(outfile)

    return arguments["o"]


includes = normalize_expr(includes, "include")
excludes = normalize_expr(excludes, "exclude", len(includes))
includes.update(excludes)

if index and not gz:
    logger.warning("Forcing envs.gz to True because envs.index is True.")
    gz = True

envs[""] = [bcftools, "filter"]
envs["_"] = infile
envs["o"] = outfile
envs["threads"] = ncores

if "O" not in envs and "output-type" not in envs and "output_type" not in envs:
    envs["O"] = "z" if gz else "v"

if keep:
    envs["soft_filter"] = "+"

if "m" not in envs and "mode" not in envs:
    envs["m"] = "+"

# bcftools can be only done once at one filter
for i, (fname, (filt, flag)) in enumerate(includes.items()):
    infile = handle_filter(infile, fname, filt, flag, i == len(includes) - 1)
