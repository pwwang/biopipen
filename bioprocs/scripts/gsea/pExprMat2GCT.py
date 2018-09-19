samples = []
ngenes  = 0
with open("{{i.expfile}}") as f:
	for line in f:
		if not samples:
			samples = line.strip().split('\t')
		elif line.strip():
			ngenes += 1

with open("{{i.expfile}}") as f, open("{{o.outfile}}", "w") as fout:
	fout.write("#1.2\n")
	fout.write("%s\t%s\n" % (ngenes, len(samples)))
	fout.write("NAME\tDescription\t%s\n" % ('\t'.join(samples)))
	f.readline() # header
	for line in f:
		line  = line.strip()
		if not line: continue
		parts = line.split('\t')
		parts.insert(1, parts[0])
		fout.write('\t'.join(parts) + '\n')
