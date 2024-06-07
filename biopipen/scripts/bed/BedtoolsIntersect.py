from pathlib import Path
from biopipen.utils.misc import run_command, dict_to_cli_args, logger

afile = {{in.afile | repr}}  # pyright: ignore # noqa: #999
bfile = {{in.bfile | repr}}  # pyright: ignore
outfile = {{out.outfile | repr}}  # pyright: ignore
envs = {{envs | repr}}  # pyright: ignore

bedtools = envs.pop("bedtools")
sort = envs.pop("sort")
chrsize = envs.pop("chrsize")
outdir = Path(outfile).parent

if chrsize and "g" in envs:
    logger.warning("Ignoring envs.g because envs.chrsize is provided.")
    envs["g"] = chrsize
elif chrsize:
    envs["g"] = chrsize

if sort:
    afile_sorted = outdir / f"{afile.stem}_sorted{afile.suffix}"
    bfile_sorted = outdir / f"{bfile.stem}_sorted{bfile.suffix}"
    run_command(f'sort -k1,1 -k2,2n -o "{afile_sorted}" "{afile}"', fg=True)
    run_command(f'sort -k1,1 -k2,2n -o "{bfile_sorted}" "{bfile}"', fg=True)
    afile = afile_sorted
    bfile = bfile_sorted

envs[""] = [bedtools, "intersect"]
envs["a"] = afile
envs["b"] = bfile
envs.setdefault("sorted", True)

if envs["sorted"] and not "g" in envs:
    raise ValueError("envs.g is required or manullay set envs.sorted to False.")

run_command(dict_to_cli_args(envs, prefix="-"), stdout=outfile)
