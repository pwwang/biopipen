colnames = []
rownames = set()
data     = {}
with open('{{i.infile}}') as f:
	for line in f:
		line  = line.strip()
		if not line: continue
		parts = line.split('\t')
		rname = parts.pop(0)
		colnames.append(rname)
		parts.pop(0) # annotation
		data[rname] = parts
		rownames = rownames | set(parts)

with open('{{o.outfile}}', 'w') as f:
	f.write("\t" + "\t".join(colnames) + "\n")
	for rname in rownames:
		tmp = ['1' if rname in data[cname] else '0' for cname in colnames]
		tmp.insert(0, rname)
		f.write("\t".join(tmp) + "\n")
