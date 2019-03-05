import re
from os import path
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs, logger
from bioprocs.utils.meme import MemeReader
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter, TsvRecord
from bioprocs.utils.parallel import Parallel

tffile   = {{i.tffile | quote}}
sfile    = {{i.sfile | quote}}
outfile  = {{o.outfile | quote}}
outdir   = {{o.outdir | quote}}
tool     = {{args.tool | quote}}
meme     = {{args.meme | quote}}
params   = {{args.params | repr}}
tfmotifs = {{args.tfmotifs | quote}}
pval     = {{args.pval | repr}}
ucsclink = {{args.ucsclink | quote}}
nthread  = {{args.nthread | repr}}

# get all motifs
logger.info('Loading motif names ...')
reader = TsvReader(tffile, cnames = False)
ncol   = len(next(reader))
reader.rewind()
if ncol == 1:
	motifns = {r[0]:r[0] for r in reader}
else:
	motifns = {r[0]:r[1] for r in reader}
logger.info('%s motif names read.', len(motifns))

# match motifs
logger.info('Matching motifs in database ...')
mnames = [m.name for m in MemeReader(tfmotifs)]
motifs = {k:v for k, v in motifns.items() if k in mnames}
logger.info('%s motifs loaded', len(motifs))

if tool == 'meme':
	cmdparams        = []
	params.thresh    = pval
	params.verbosity = 4
	for motif, name in motifs.items():
		params.oc    = path.join(outdir, name + '.' + re.sub(r'[^\w_]', '', motif))
		params.motif = motif
		params._     = [tfmotifs, sfile]
		cmdparams.append((meme, cmdargs(params, dash = '--', equal = ' ')))
	Parallel(nthread, raiseExc = True).run('{} {}', cmdparams)

	writer = TsvWriter(outfile)
	writer.cnames = [
		"CHR", "START", "END", "NAME", "SCORE", "STRAND", "MOTIF", "SEQ", "STARTONSEQ",
		"STOPONSEQ", "RAWSCORE", "PVAL", "QVAL", "MATCHEDSEQ", "UCSCLINK"
	]
	writer.writeHead(callback = lambda cnames: "#" + "\t".join(cnames))

	def rowfactory(r):
		if 'p-value' not in r:
			return None
		r.PVAL       = float(r['p-value'])
		if r.PVAL >= pval:
			return None
		r.RAWSCORE = r.score
		try:
			r.SCORE = int(float(r.score) * 10)
		except TypeError:
			r.SCORE = 0
		r.STRAND   = r.strand
		r.MOTIF    = r.motif_id
		# split motif_alt_id
		# GENE or GENE::chr1:111-222 or ::chr1:111-222 or chr1:111-222
		r.SEQ = r.sequence_name
		if '::' in r.SEQ:
			seqname, seqcoord = r.SEQ.split('::', 1)
			r.CHR, start, end = re.split(r'[:-]', seqcoord)
			start = int(start)
			end   = int(end)
		elif ':' in r.SEQ:
			seqname, seqcoord = r.SEQ, r.SEQ
			r.CHR, start, end = re.split(r'[:-]', seqcoord)
			start = int(start)
			end   = int(end)
		else:
			seqname, seqcoord = r.SEQ, ''
			r.CHR, start, end = '-', 0, 0
		r.NAME       = seqname
		r.STARTONSEQ = r.start
		r.STOPONSEQ  = r.stop
		r.START      = start + int(r.start)
		r.END        = start + int(r.stop)
		r.QVAL       = r['q-value']
		r.MATCHEDSEQ = r.matched_sequence
		r.UCSCLINK   = ucsclink.format('{chrom}:{start}:{end}'.format(
			chrom = r.CHR,
			start = start,
			end   = end
		))
		return r

	for motif, name in motifs.items():
		retfile = path.join(outdir, name + '.' + re.sub(r'[^\w_]', '', motif), 'fimo.tsv')
		# motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence
		logger.info('- Merging %s ...', retfile)
		reader = TsvReader(
			retfile, 
			#cnames  = lambda header: header[1:].strip().split('\t'),
			cnames  = True,
			row     = rowfactory,
			comment = None)
		for r in reader:
			if not r:
				continue
			r.NAME   = name + '::' + r.NAME
			writer.write(r)
	writer.close()
else:
	raise ValueError('Unsupported tool: {}'.format(tool))




