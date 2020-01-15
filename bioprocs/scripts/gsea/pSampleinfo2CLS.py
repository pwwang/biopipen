
from bioprocs.utils.sampleinfo import SampleInfo

infile    = {{i.sifile | quote}}
outfile   = {{o.outfile | quote}}
saminfo   = SampleInfo(infile)
groups    = saminfo.select(get = 'Group')
unigroups = list(set(groups))

with open(outfile, "w") as f:
	f.write("%s\t%s\t1\n" % (saminfo.nrow, len(unigroups)))
	f.write("# %s\n" % (' '.join(unigroups)))
	f.write(' '.join(groups) + '\n')
