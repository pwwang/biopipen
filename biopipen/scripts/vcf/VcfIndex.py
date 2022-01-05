from os import path

import cmdy
from biopipen.utils.reference import tabix_index

infile = {{in.infile | repr}}
outfile = {{out.outfile | repr}}
outidx = {{out.outidx | repr}}
tabix = {{envs.tabix | repr}}
ncores = {{envs.ncores | repr}}

outfile_with_index = tabix_index(infile, "vcf", outfile.parent, tabix)
if path.samefile(infile, outfile_with_index):
    cmdy.ln(s=infile, _=outfile)
    cmdy.ln(s=infile + ".tbi", _=outidx)
