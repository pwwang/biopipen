from pyppl import Box
from bioprocs.utils.shell import grep, zcat

infile  = {{ i.infile | quote}}
outfile = {{ o.outfile | quote}}
params  = {{ args.params | repr}}
keyword = {{ args.keyword | repr}}
if not keyword:
	raise ValueError('A keyword (args.keyword) is required.')

params._stdout = outfile
if infile.endswith('.gz'):
	zcat(infile).pipe().grep(_ = keyword, **params)
else:
	grep(_ = [keyword, infile], **params)