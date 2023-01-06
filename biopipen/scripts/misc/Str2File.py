instr = {{in.str | repr}}  # pyright: ignore
name = {{repr(in.name or envs.name)}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore

with open(outfile, "wt") as fout:
    fout.write(instr)
