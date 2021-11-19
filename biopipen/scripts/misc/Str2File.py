instr = {{in.str | repr}}
name = {{repr(in.name or envs.name)}}
outfile = {{out.outfile | quote}}

with open(outfile, "wt") as fout:
    fout.write(instr)
