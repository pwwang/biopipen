from diot import Diot
from bioprocs.utils.shell import grep, zcat

infile  = {{ i.infile | quote}}
outfile = {{ o.outfile | quote}}
params  = {{ args.params | repr}}
keyword = {{ args.keyword | repr}}
if not keyword:
	raise ValueError('A keyword (args.keyword) is required.')

# params._stdout = outfile
if infile.endswith('.gz'):
	zcat(infile).p | grep(_ = keyword, **params) ^ outfile
else:
	grep(_ = [keyword, infile], **params)
