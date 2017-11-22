{{txtFilter}}

txtFilter({{in.infile | quote}}, {{out.outfile | quote}}, cols = [0, 2, "V4"], rfilter=lambda row: float(row[1]) > 2 and float(row[2]) > 2, header = True, skip = 2, delimit = '|')