instr: str = {{in.str | quote}}  # pyright: ignore  # noqa
name = {{repr(in.name or envs.name)}}  # pyright: ignore
outfile: str = {{out.outfile | quote}}  # pyright: ignore

with open(outfile, "wt") as fout:
    fout.write(instr)
