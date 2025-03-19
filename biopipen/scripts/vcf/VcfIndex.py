from pathlib import Path
from os import path

from biopipen.utils.reference import tabix_index
from biopipen.utils.misc import run_command

infile: str = {{in.infile | quote}}  # pyright: ignore  # noqa
outfile = Path({{out.outfile | quote}})  # pyright: ignore
outidx = {{out.outidx | repr}}  # pyright: ignore
tabix: str = {{envs.tabix | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore

outfile_with_index = tabix_index(infile, "vcf", outfile.parent, tabix)
if path.samefile(infile, outfile_with_index):
    run_command(["ln", "-s", infile, outfile], fg=True)
    run_command(["ln", "-s", infile + ".tbi", outidx], fg=True)
