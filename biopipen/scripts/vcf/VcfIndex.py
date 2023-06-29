from pathlib import Path
from os import path

from biopipen.utils.reference import tabix_index
from biopipen.utils.misc import run_command

infile = {{in.infile | repr}}  # pyright: ignore
outfile = Path({{out.outfile | repr}})  # pyright: ignore
outidx = {{out.outidx | repr}}  # pyright: ignore
tabix = {{envs.tabix | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore

outfile_with_index = tabix_index(infile, "vcf", outfile.parent, tabix)
if path.samefile(infile, outfile_with_index):
    run_command(["ln", "-s", infile, outfile], fg=True)
    run_command(["ln", "-s", infile + ".tbi", outidx], fg=True)
