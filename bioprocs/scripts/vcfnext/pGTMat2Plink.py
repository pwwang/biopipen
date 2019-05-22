from os import path
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter
from bioprocs.utils import logger, shell2 as shell

infile   = {{i.infile | quote}}
metafile = {{i.metafile | quote}}
outdir   = {{o.outdir | quote}}
plink    = {{args.plink | quote}}
keeptxt  = {{args.keeptxt | repr}}
chrmaps  = {{args.chrmaps | str}}
prefix   = path.join(outdir, {{i.infile | fn2 | quote}})
tpedfile  = prefix + ".tped"
tfamfile  = prefix + ".tfam"

shell.load_config(plink = plink)

# column names could be:
# FID, IID, PID, MID, Sex, Pheno
if metafile:
	logger.info('Reading metafile ...')
	metadata = dict(TsvReader(
		metafile,
		cnames = True,
		row    = lambda r: tuple((r.IID, r))
	).dump())
else:
	metadata = None

logger.info('Reading genotype matrix ...')
# snp1 gt1s1 gt1s2 ...
inreader = TsvReader(infile, cnames = True)
samples  = inreader.meta[1:]

logger.info('Writing tfam file ...')
tfamWriter = TsvWriter(tfamfile)
tfamWriter.meta = ['FID', 'IID', 'PID', 'MID', 'Sex', 'Pheno']
#tfamWriter.writeHead(callback = lambda meta: '#' + '\t'.join(meta))
uniqueIDs = {}
if not metadata:
	for s in samples:
		tfamWriter.write([s, s, '0', '0', 'other', '-9'])
else:
	for s in samples:
		fid     = metadata[s].FID if s in metadata and 'FID' in metadata[s] else s
		idcheck = fid + ' ' + s
		if idcheck in uniqueIDs:
			raise ValueError('Duplicated ID {!r}'.format(idcheck))
		uniqueIDs[idcheck] = 1
		tfamWriter.write([
			fid,
			s,
			(metadata[s].PID or '0') if s in metadata and 'PID' in metadata[s] else '0',
			(metadata[s].MID or '0') if s in metadata and 'MID' in metadata[s] else '0',
			(metadata[s].Sex or 'other') if s in metadata and 'Sex' in metadata[s] else 'other',
			(metadata[s].Pheno or '-9') if s in metadata and 'Pheno' in metadata[s] else '-9'
		])
tfamWriter.close()

def getCompondGT(gt, ref, alt):
	compGTs = {
		"0": ref + ' ' + ref,
		"1": ref + ' ' + alt,
		"2": alt + ' ' + alt
	}
	return compGTs.get(gt, '0 0')

logger.info('Writing tped file ...')
tpedWriter = TsvWriter(tpedfile)
for r in inreader:
	(chrom, pos, _, ref, alt) = r[0].split('_')
	if chrom.startswith('chr'):
		chrom = chrom[3:]
	chrom = chrmaps.get(chrom, chrom)
	tpedWriter.write([chrom, r[0], 0, pos] + [getCompondGT(gt, ref, alt) for gt in r.values()[:]])
tpedWriter.close()

logger.info("Converting using plink ...")
shell.fg.plink({'make-bed': True}, tfile = prefix, out = prefix)

if not keeptxt:
	shell.rm_rf(tpedfile)
	shell.rm_rf(tfamfile)
