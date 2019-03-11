from collections import OrderedDict
from bioprocs.utils import shell, FileConn
from bioprocs.utils.reference import vcfIndex
from bioprocs.utils.parallel import Parallel

{% python from bioprocs.utils import alwaysList %}
infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
sortby   = {{args.sortby | quote}}
picard   = {{args.picard | quote}}
tabix    = {{args.tabix | quote}}
nthread  = {{args.nthread | repr}}
gsize    = {{args.gsize | quote}}
chrorder = {{args.chrorder | alwaysList | repr}}
tool     = {{args.tool | quote}}

def getContigs(infile):
	ret = OrderedDict()
	with FileConn(infile) as fh:
		for line in fh:
			if not line.startswith('#'):
				break
			##contig=<ID=chr1,length=249250621>
			if not line.startswith('##contig='):
				break
			items = line.rstrip('>\n')[10:].split(',')
			tmp   = {}
			for item in items:
				key, val = item.split('=', 1)
				tmp[key] = val
				if key == 'ID':
					ret[val] = tmp
	return ret

def splitChroms(infile, chroms):
	fchrs = {}
	for chrom in chroms:
		fchrs[chrom] = open(outfile + '.' + chrom + '.split', 'w')
	with FileConn(infile) as fh:
		for line in fh:
			if line.startswith('#'):
				continue
			chrom = line.split('\t', 1)[0]
			if not chrom in fchrs:
				continue
			fchrs[chrom].write(line)
	for chrom in chroms:
		fchrs[chrom].close()

def cleanChroms(para = None):
	from glob import glob
	para = para or Parallel(nthread)
	para.run(shell.rmrf, [(fn, ) for fn in glob(outfile + '.*.s*')]) # .split, .sorted
	del para

def copyHeader(hasContigs, contigs = None):
	with FileConn(infile) as fh, open(outfile, 'w') as fout:
		for line in fh:
			if not line.startswith('#'):
				break
			if line.startswith('##contig='):
				if hasContigs and not contigs:
					fout.write(line)
				else:
					continue
			elif line.startswith('##'):
				fout.write(line)
			elif not contigs: # startswith '#'
				fout.write(line)
			else:
				for chrom, contig in contigs.items():
					fout.write('##contig=<ID={},length={}>\n'.format(contig['ID'], contig['length']))
				fout.write(line)

def run_sort_nochrorder_nocontigs(sortby):
	copyHeader(hasContigs = False, contigs = None)
	if sortby == 'coord':
		params = dict(k = ['1,1', '2,2n'], _stdout_ = outfile)
	else:
		params = dict(k = ['3,3', '1,1', '2,2n'], _stdout_ = outfile)
	shell.acat_(infile).pipe().grep(v = '^#').pipe(duplistkey = True).sort(**params).run()

def run_sort_contigs(sortby, hasContigs, contigs): # no chrorder
	copyHeader(hasContigs = hasContigs, contigs = contigs if chrorder or not hasContigs else None)
	splitChroms(infile, contigs.keys())
	shsort = shell.Shell(duplistkey = True).sort
	if sortby == 'coord':
		sortchr = lambda chrom: shsort(k = '2,2n', _ = outfile + '.' + chrom + '.split', _stdout = outfile + '.' + chrom + '.sorted').run()
	else:
		sortchr = lambda chrom: shsort(k = ['3,3', '2,2n'], _ = outfile + '.' + chrom + '.split', _stdout = outfile + '.' + chrom + '.sorted').run()
	para = Parallel(nthread)
	para.run(sortchr, [(chrom, ) for chrom in contigs.keys()])
	for chrom in contigs.keys():
		ofile = outfile + '.' + chrom + '.sorted'
		with open(ofile) as fin, open(outfile, 'a') as fout:
			fout.write(fin.read())
	cleanChroms(para)

def chrorder_to_contigs(chrorder, gsize):
	contigs = OrderedDict()
	tmpcontigs = {}
	with open(gsize) as f:
		for line in f:
			ID, length = line.strip().split('\t')[:2]
			tmpcontigs[ID] = dict(ID = ID, length = length)
	for chrom in chrorder:
		if not chrom in tmpcontigs:
			continue
		contigs[chrom] = tmpcontigs[chrom]
	return contigs

def run_sort():
	contigs = getContigs(infile)
	hasContigs = bool(contigs)
	if not chrorder and not contigs:
		run_sort_nochrorder_nocontigs(sortby = sortby)
	else:
		if chrorder:
			contigs = chrorder_to_contigs(chrorder, gsize)
		run_sort_contigs(sortby = sortby, hasContigs = hasContigs, contigs = contigs)

def run_picard():
	if sortby != 'coord':
		raise ValueError('Picard does not sort by name, use sort instead.')
	shpicard = shell.Shell({'picard': picard}, subcmd = True, dash = '', equal = '=').picard
	# compose a dict file from chrorder
	if chrorder:
		contigs = chrorder_to_contigs(chrorder, gsize)
		dictfile = outfile + '.dict'
		with open(dictfile, 'w') as f:
			f.write("@HD\tVN:1.5\tSO:unsorted\n")
			f.write("".join(["@SQ\tSN:{}\tLN:{}\n".format(chrom, contig['length']) for chrom, contig in contigs.items()]))
		shpicard.SortVcf(INPUT = infile, OUTPUT = outfile, SEQUENCE_DICTIONARY = dictfile).run()
	else:
		shpicard.SortVcf(INPUT = infile, OUTPUT = outfile).run()

tools = dict(
	sort   = run_sort,
	picard = run_picard
)

try:
	tools[tool]()
except KeyError:
	raise ValueError('Tool {!r} not supported yet.'.format(tool))
except:
	raise
