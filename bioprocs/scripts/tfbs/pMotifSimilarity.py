from os import path
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter
from bioprocs.utils.parallel import Parallel, distribute
from bioprocs.utils.meme import MemeReader, MemeWriter

mfile1  = {{ i.mfile1 | quote}}
mfile2  = {{ i.mfile2 | quote}}
outfile = {{ o.outfile | quote}}
outdir  = {{ o.outdir | quote}}
qval    = {{ args.qval | repr}}
tomtom  = {{ args.tomtom | quote}}
params  = {{ args.params | repr}}
nthread = {{ args.nthread | repr }}
mfile2  = mfile2 or mfile1

if nthread > 1:
	reader  = MemeReader(mfile1)
	count   = sum(1 for r in reader)
	joblist = distribute(count, nthread)
	cmdps   = []
	ocdirs  = []
	reader.rewind()
	for i in range(nthread):
		qfile = path.join(outdir, 'query-{}.meme'.format(i+1))
		ocdir = path.join(outdir, 'query-{}.tomtom'.format(i+1))
		ocdirs.append(ocdir)
		writer = MemeWriter(qfile)
		writer.meta = reader.meta
		writer.writeMeta()
		for _ in range(joblist[i]):
			try:
				writer.write(reader.next())
			except StopIteration:
				break
		writer.close()
		thparams = params.copy()
		thparams[""] = [qfile, mfile2]
		thparams.thresh = qval
		thparams.oc  = ocdir
		cmdps.append((tomtom, cmdargs(thparams, dash = '-', equal = ' ')))
	reader.close()
	Parallel(nthread, raiseExc = True).run('{} {}', cmdps)

	writer = TsvWriter(outfile)
	reader = TsvReader(path.join(ocdirs[0], 'tomtom.txt'), comment = '##', cnames = lambda header: header[1:].strip().split("\t"))
	writer.cnames = reader.cnames
	writer.writeHead(lambda cnames: "#" + "\t".join(cnames))
	reader.close()
	for ocdir in ocdirs:
		reader = TsvReader(path.join(ocdir, 'tomtom.txt'), comment = '##', cnames = lambda header: header[1:].strip().split("\t"))
		for r in reader:
			writer.write(r)
		reader.close()
	writer.close()
else:
	params[""] = [mfile1, mfile2]
	params.thresh = qval
	params.oc = outdir

	cmd = '{tomtom} {params}'.format(
		tomtom = tomtom,
		params = cmdargs(params, dash = '-', equal = ' ')
	)
	runcmd(cmd)
