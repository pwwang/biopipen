{{ txtSampleinfo }}
info, cols = txtSampleinfo("{{in.sifile}}")
groups = list(set(cols['Group']))

with open("{{out.outfile}}", "w") as f:
	f.write("%s\t%s\t1\n" % (len(info), len(groups)))
	f.write("# %s\n" % (' '.join(groups)))
	clss = []
	for sample, sinfo in info.items():
		clss.append(sinfo['Group'])
	f.write(' '.join(clss) + '\n')