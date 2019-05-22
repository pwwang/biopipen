from os import path, makedirs, remove
from shutil import rmtree, copyfile, move
from sys import stderr
from pyppl import Box
from bioprocs.utils import log2pyppl
from bioprocs.utils.shell2 import bgzip
import vcf

infile  = {{i.infile | quote}}
filters = {{args.filters | repr}}
fname   = {{args.filters | quote}}
gz      = {{args.gz | repr}}
outfile = {{o.outfile | quote}}
keep    = {{args.keep | repr}}
if gz:
	outfile = outfile[:-3]

"""
Builtin filters
"""
builtin_filters = Box(
	SNPONLY   = lambda r, s, q = None: len(r.REF) != 1 or any(len(a) != 1 for a in r.ALT),
	BIALTONLY = lambda r, s, q = None: len(r.ALT) != 1,
	QUAL      = lambda r, s, q: not r.QUAL or r.QUAL < q,
	AF        = lambda r, s, q: r.INFO['AF'][0] < q
)
builtin_descs = Box(
	SNPONLY   = lambda x: 'keep SNPs only',
	BIALTONLY = lambda x: 'keep bi-allelic mutations only',
	QUAL      = lambda x: 'keep sites with QUAL >= %s' % x,
	AF        = lambda x: 'keep sites with AF >= %s' % x
)
desc_prefix = 'Created by bioprocs.vcf.pVcfFilter: '

filters     = filters or {}
realfilters = {}
descs       = {}
for fname, ffunc in filters.items():
	if fname.startswith('!') and fname[1:] in builtin_filters:
		key = 'NOT' + fname[1:] + ['', str(ffunc)][int(bool(ffunc) or 0)]
		realfilters[key] = lambda r, s, fname = fname, ffunc = ffunc: \
			not builtin_filters[fname](r, s, ffunc) if ffunc not in [None, False] \
				else not builtin_filters[fname](r, s)
		descs[key] = desc_prefix + builtin_descs[fname[1:]](ffunc)
	elif fname in builtin_filters:
		key = [fname, fname + str(ffunc)][int(bool(ffunc) or 0)]
		realfilters[key] = lambda r, s, fname = fname, ffunc = ffunc: \
			builtin_filters[fname](r, s, ffunc) if ffunc not in [None, False] \
				else builtin_filters[fname](r, s)
		descs[key] = desc_prefix + builtin_descs[fname](ffunc)
	else:
		realfilters[fname] = ffunc if callable(ffunc) else eval(ffunc)
		descs[fname] = desc_prefix + fname

reader = vcf.Reader(filename=infile)
for fname, fdesc in descs.items():
	reader.filters[fname] = vcf.parser._Filter(id = fname, desc = fdesc)
writer = vcf.Writer(open(outfile, 'w'), reader)

while True:
	nerror = 0
	try:
		record = next(reader)
		for fname, ffunc in realfilters.items():
			if ffunc(record, record.samples):
				record.FILTER = record.FILTER or []
				record.FILTER.append(fname)
		if keep or not record.FILTER:
			writer.write_record(record)
	except StopIteration:
		break
	except Exception:
		nerror += 1
		import traceback
		tb = traceback.format_exc()
		log2pyppl(str(tb), 'ERROR')
		if nerror < 100:
			continue
		raise

writer.close()

if gz:
	bgzip(outfile)


