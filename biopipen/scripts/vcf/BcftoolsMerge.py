from biopipen.utils.reference import tabix_index
from biopipen.utils.misc import logger
from biopipen.scripts.vcf.bcftools_utils import run_bcftools

infiles: list = {{in.infiles | each: as_path}}  # pyright: ignore # noqa: E999
outfile = {{out.outfile | repr}}  # pyright: ignore
joboutdir = {{job.outdir | repr}}  # pyright: ignore
envs: dict = {{envs | dict | repr}}  # pyright: ignore

bcftools = envs.pop("bcftools")
tabix = envs.pop("tabix")
ncores = envs.pop("ncores")
gz = envs.pop("gz")
index = envs.pop("index")

envs.setdefault("force-single", True)
envs.setdefault("missing-to-ref", True)

if index and not gz:
    logger.warning("Forcing envs.gz to True because envs.index is True.")
    gz = True

if "O" not in envs and "output-type" not in envs and "output_type" not in envs:
    envs["O"] = "z" if gz else "v"

envs[""] = [bcftools, "merge"]
envs["o"] = outfile
envs["threads"] = ncores
envs["_"] = infiles

run_bcftools(envs, bcftools=bcftools, index=index, tabix=tabix)
