from os import path

import cmdy
from biopipen.utils.reference import tabix_index

infile = {{in.infile | repr}}
annfile = {{(in.annfile or envs.annfile) | repr}}
outfile = {{out.outfile | repr}}
joboutdir = {{job.outdir | repr}}
bcftools = {{envs.bcftools | repr}}
tabix = {{envs.tabix | repr}}
ncores = {{envs.ncores | repr}}
cols = {{envs.cols | repr}}
header = {{envs.header | repr}}
args = {{envs.args | repr}}

args["_exe"] = bcftools
args["_"] = tabix_index(infile, "vcf")
args["o"] = outfile
args["threads"] = ncores

if annfile:
    abname = path.basename(annfile)
    ext = path.splitext(
        abname[:-3] if abname.endswith('.gz') else abname
    )[-1][1:]
    args["a"] = tabix_index(annfile, ext, tabix)

if cols and isinstance(cols, list):
    args["c"] = ",".join(cols)

if header:
    if not isinstance(header, list):
        header = [header]

    headerfile = path.join(joboutdir, "header.txt")
    with open(headerfile, "w") as fh:
        for head in header:
            fh.write(f"{head}\n")
    args["h"] = headerfile

cmd = cmdy.bcftools.annotate(**args).h()
print("Running:")
print("-------")
print(cmd.strcmd)
cmd.fg().run()
