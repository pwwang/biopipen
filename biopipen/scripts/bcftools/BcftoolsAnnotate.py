from os import path

from biopipen.utils.reference import tabix_index
from biopipen.utils.misc import dict_to_cli_args, run_command

infile = {{in.infile | repr}}  # pyright: ignore
annfile = {{(in.annfile or envs.annfile) | repr}}  # pyright: ignore
outfile = {{out.outfile | repr}}  # pyright: ignore
joboutdir = {{job.outdir | repr}}  # pyright: ignore
bcftools = {{envs.bcftools | repr}}  # pyright: ignore
tabix = {{envs.tabix | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore
cols = {{envs.cols | repr}}  # pyright: ignore
header = {{envs.header | repr}}  # pyright: ignore
args = {{envs.args | repr}}  # pyright: ignore

args[""] = bcftools
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

run_command(dict_to_cli_args(args, dashify=True), fg=True)
