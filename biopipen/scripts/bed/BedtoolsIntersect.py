from pathlib import Path
from biopipen.utils.misc import run_command, dict_to_cli_args, logger

afile = Path({{in.afile | quote}})  # pyright: ignore # noqa: #999
bfile = Path({{in.bfile | quote}})  # pyright: ignore
outfile: str = {{out.outfile | quote}}  # pyright: ignore
envs: dict = {{envs | dict}}  # pyright: ignore

bedtools = envs.pop("bedtools")
sort = envs.pop("sort")
chrsize = envs.pop("chrsize")
postcmd = envs.pop("postcmd", None)
outdir = Path(outfile).parent

if chrsize and "g" in envs:
    logger.warning("Ignoring envs.g because envs.chrsize is provided.")
    envs["g"] = Path(chrsize).expanduser()
elif chrsize:
    envs["g"] = Path(chrsize).expanduser()

if sort:
    afile_sorted = outdir / f"{afile.stem}_sorted{afile.suffix}"
    bfile_sorted = outdir / f"{bfile.stem}_sorted{bfile.suffix}"
    run_command(
        [bedtools, "sort", "-g", envs["g"], "-i", afile],
        stdout=afile_sorted,
    )
    run_command(
        [bedtools, "sort", "-g", envs["g"], "-i", bfile],
        stdout=bfile_sorted,
    )
    afile = afile_sorted
    bfile = bfile_sorted

envs[""] = [bedtools, "intersect"]
envs["a"] = afile
envs["b"] = bfile
envs.setdefault("sorted", True)

if envs["sorted"] and not "g" in envs:
    raise ValueError("envs.g is required or manullay set envs.sorted to False.")

if postcmd:
    ofile = Path(outfile).with_suffix(".prior.bt")
    run_command(dict_to_cli_args(envs, prefix="-"), stdout=ofile)
    postcmd_file = outdir / "_postcmd.sh"
    postcmd_file.write_text(postcmd)
    run_command(
        ["bash", postcmd_file],
        env={"infile": ofile, "outfile": outfile, "outdir": outdir},
        fg=True,
    )
else:
    run_command(dict_to_cli_args(envs, prefix="-"), stdout=outfile)
