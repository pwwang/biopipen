import cmdy

infile1 = {{in.infile1 | repr}}  # pyright: ignore
infile2 = {{in.infile2 | repr}}  # pyright: ignore
outfile = {{out.outfile | repr}}  # pyright: ignore
bcftools = {{envs.bcftools | repr}}  # pyright: ignore
gz = {{envs.gz | repr}}  # pyright: ignore
index = {{envs.index | repr}}  # pyright: ignore
collapse = {{envs.collapse | repr}}  # pyright: ignore

if index:
    gz = True

args = {
    "c": collapse,
    "O": "z" if gz else "v",
    "o": outfile,
    "write": 1,
    "_": ["-n=2", infile1, infile2],
    "_exe": bcftools,
}

cmd = cmdy.bcftools.isec(**args).hold()
print("  running:")
print("  ", cmd.strcmd)
cmd.run(wait=True)

if index:
    cmd = cmdy.bcftools.index(_=outfile, t=True, _exe=bcftools).hold()
    print("  running:")
    print("  ", cmd.strcmd)
    cmd.run(wait=True)
