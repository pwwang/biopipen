{{txtTransform}}

txtTransform({{in.infile | quote}}, {{out.outfile | quote}}, cols = [0, 2, "V4"], transform=lambda row: [str(float(r) + 1) if i == 1 else r for i, r in enumerate(row)], header = True, skip = 2, delimit = '|')