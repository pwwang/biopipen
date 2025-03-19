
infile: str = {{in.infile | quote}}  # pyright: ignore  # noqa
outfile: str = {{out.outfile | quote}}  # pyright: ignore

with open(infile) as fin, open(outfile, "w") as fout:
    for line in fin:
        if line.startswith("#") or line.startswith("Hugo_Symbol"):
            fout.write(line)
        else:
            cols = line.split("\t")
            if not cols[4].startswith("chr"):
                cols[4] = f"chr{cols[4]}"
            # "\n" at the last col kept
            fout.write("\t".join(cols))
