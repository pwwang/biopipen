import cmdy

infile = {{in.infile | quote}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore
bcftools = {{envs.bcftools | quote}}  # pyright: ignore
gz = {{envs.gz | repr}}  # pyright: ignore
args = {{envs.args | repr}}  # pyright: ignore
tmpdir = {{envs.tmpdir | quote}}  # pyright: ignore
index = {{envs.index | repr}}  # pyright: ignore

args["_exe"] = bcftools
args["_"] = infile
args["o"] = outfile
args["O"] = "z" if gz or index else "v"

cmdy.bcftools.sort(**args).fg()

if index:
    cmdy.bcftools.index(outfile).fg()
