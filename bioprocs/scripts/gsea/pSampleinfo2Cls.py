{{ txtSampleinfo }}
saminfo = txtSampleinfo("{{in.sifile}}")
groups = [g[0] for g in saminfo[["", "Group"]].data]

with open("{{out.outfile}}", "w") as f:
	f.write("%s\t%s\t1\n" % (saminfo.nrow, len(set(groups))))
	f.write("# %s\n" % (' '.join(list(set(groups)))))
	
	f.write(' '.join(groups) + '\n')