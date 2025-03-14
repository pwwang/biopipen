from os import path
from contextlib import suppress
from pathlib import PosixPath  # noqa: F401

from biopipen.utils.reference import tabix_index
from biopipen.utils.misc import logger
from biopipen.scripts.vcf.bcftools_utils import run_bcftools

infile: str = {{in.infile | quote}}  # pyright: ignore # noqa: E999
annfile: str = {{in.annfile | quote}}  # pyright: ignore
outfile: str = {{out.outfile | quote}}  # pyright: ignore
joboutdir: str = {{job.outdir | quote}}  # pyright: ignore
envs: dict = {{envs | dict | repr}}  # pyright: ignore

bcftools = envs.pop("bcftools")
tabix = envs.pop("tabix")
ncores = envs.pop("ncores")
columns = envs.pop("columns")
remove = envs.pop("remove")
header = envs.pop("header")
gz = envs.pop("gz")
index = envs.pop("index")

if isinstance(columns, list):
    columns = ",".join(columns)

if "c" in envs:
    logger.warning(r"Ignoring envs\[c], use envs\[columns] instead.")
    del envs["c"]

if isinstance(remove, list):
    remove = ",".join(remove)

if "x" in envs:
    logger.warning(r"Ignoring envs\[x], use envs\[remove] instead.")
    del envs["x"]

envs_has_annfile = "a" in envs or "annotations" in envs
headerfile = path.join(joboutdir, "header.txt")
if header:
    with open(headerfile, "w") as fh:
        fh.writelines(header)

if annfile and envs_has_annfile:
    logger.warning(
        r"Ignoring envs\[a/annotations] because in.annfile is provided."
    )
    with suppress(KeyError):
        del envs["a"]
    with suppress(KeyError):
        del envs["annotations"]
elif not annfile and envs_has_annfile:
    annfile = envs.pop("annotations", None) or envs.pop("a", None)


if index and not gz:
    logger.warning("Forcing envs.gz to True because envs.index is True.")
    gz = True

envs[""] = [bcftools, "annotate"]
envs["o"] = outfile
envs["threads"] = ncores

if "O" not in envs and "output-type" not in envs and "output_type" not in envs:
    envs["O"] = "z" if gz else "v"

if columns:
    envs["columns"] = columns
    if not annfile:
        raise ValueError(
            "envs.columns specified but no in.annfile/envs.annfile provided."
        )
    envs["_"] = tabix_index(infile, "vcf", tabix=tabix)

if remove:
    envs["remove"] = remove
    # no need to index it
    envs["_"] = infile

if "columns" not in envs and "remove" not in envs:
    logger.warning(
        "No columns/remove specified, no columns will be carried over or removed."
    )

if annfile:
    envs["annotations"] = tabix_index(annfile, "vcf", tabix=tabix)

if header:
    envs["header_lines"] = headerfile

run_bcftools(envs, bcftools=bcftools, index=index, tabix=tabix)
