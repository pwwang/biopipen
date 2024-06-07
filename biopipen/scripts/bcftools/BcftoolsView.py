from contextlib import suppress
# In case there are paths passed to envs
from pathlib import PosixPath  # noqa: F401
from shutil import copy2

from biopipen.utils.misc import run_command, dict_to_cli_args, logger
from biopipen.utils.reference import tabix_index

infile = {{in.infile | repr}}  # pyright: ignore # noqa: #999
regions_file = {{in.regions_file | repr}}  # pyright: ignore
samples_file = {{in.samples_file | repr}}  # pyright: ignore
outfile = {{out.outfile | repr}}  # pyright: ignore
envs: dict = {{envs | dict | repr}}  # pyright: ignore

bcftools = envs.pop("bcftools")
tabix = envs.pop("tabix")
ncores = envs.pop("ncores")
gz = envs.pop("gz")
index = envs.pop("index")

if regions_file:
    if "R" in envs or "regions_file" in envs or "regions-file" in envs:
        logger.warning(
            "Ignoring envs\[regions_file/regions-file/R] "
            "because in.regionsfile is provided."
        )
        with suppress(KeyError):
            del envs["regions_file"]
        with suppress(KeyError):
            del envs["regions-file"]
        with suppress(KeyError):
            del envs["R"]
elif "R" in envs or "regions_file" in envs or "regions-file" in envs:
    regions_file = (
        envs.pop("regions_file", None)
        or envs.pop("regions-file", None)
        or envs.pop("R", None)
    )

if samples_file:
    if "S" in envs or "samples_file" in envs or "samples-file" in envs:
        logger.warning(
            "Ignoring envs[samples_file/samples-file/S] "
            "because in.samples_file is provided."
        )
        with suppress(KeyError):
            del envs["samples_file"]
        with suppress(KeyError):
            del envs["samples-file"]
        with suppress(KeyError):
            del envs["S"]
elif "S" in envs or "samples_file" in envs or "samples-file" in envs:
    samples_file = (
        envs.pop("samples_file", None)
        or envs.pop("samples-file", None)
        or envs.pop("S", None)
    )

if index and not gz:
    logger.warning("Forcing envs.gz to True because envs.index is True.")
    gz = True

if "O" not in envs and "output-type" not in envs and "output_type" not in envs:
    envs["O"] = "z" if gz else "v"

envs[""] = [bcftools, "view"]
envs["_"] = tabix_index(infile, "vcf", tabix=tabix)
envs["o"] = outfile
envs["threads"] = ncores
envs["regions_file"] = regions_file
envs["samples_file"] = samples_file

if index:
    # requires bcftools 1.20+
    # '--write-index tbi' not working
    # it has to be '--write-index=tbi'
    envs["write_index=tbi"] = True

try:
    run_command(dict_to_cli_args(envs, dashify=True), fg=True)
finally:
    run_command(["which", bcftools], fg=True)
    run_command([bcftools, "version"], fg=True)
