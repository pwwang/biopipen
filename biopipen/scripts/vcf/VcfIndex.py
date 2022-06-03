from pathlib import Path
from os import path

import cmdy
from biopipen.utils.reference import tabix_index

infile = {{in.infile | repr}}  # pyright: ignore
outfile = Path({{out.outfile | repr}})  # pyright: ignore
outidx = {{out.outidx | repr}}  # pyright: ignore
tabix = {{envs.tabix | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore

outfile_with_index = tabix_index(infile, "vcf", outfile.parent, tabix)
if path.samefile(infile, outfile_with_index):
    cmdy.ln(s=infile, _=outfile)
    cmdy.ln(s=infile + ".tbi", _=outidx)
