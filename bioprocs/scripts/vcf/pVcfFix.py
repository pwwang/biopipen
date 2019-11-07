import gzip
from os import path, environ
from cyvcf2 import VCF, Writer
from pyppl import Box
from difflib import SequenceMatcher
from bioprocs.utils import shell2 as shell

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
errfile = {{job.errfile | quote}}
# currently supported:
# clinvarLink  = True,
# addChr       = True,
# headerInfo   = {id1: {Number:1, Type:String, Description: xxx}},
# headerContig = {ref: xxx.fai/xxx.dict, notfound: drop}
# headerFormat = {id1: {Number:1, Type:String, Description: xxx}},
# headerFilter = {filter1: descirption1, filter2: desc2, ...}

fixes   = {{args.fixes | repr}}
nthread = {{args.nthread | repr}}
ref     = {{args.ref | quote}}
reffai  = ref + '.fai'
refdict = ref[:-3] + '.dict'

environ['OPENBLAS_NUM_THREADS'] = str(nthread)
environ['OMP_NUM_THREADS']      = str(nthread)
environ['NUMEXPR_NUM_THREADS']  = str(nthread)
environ['MKL_NUM_THREADS']      = str(nthread)
import numpy

inverse            = fixes.get('_inverse', False)
fixes.clinvarLink  = fixes.get('clinvarLink', False)
fixes.addChr       = fixes.get('addChr', False)
fixes.addAF        = fixes.get('addAF', False)
fixes.tumorpos     = fixes.get('tumorpos', False)
fixes.headerInfo   = fixes.get('headerInfo', False)
fixes.headerContig = fixes.get('headerContig', False)
fixes.headerFormat = fixes.get('headerFormat', False)
fixes.headerFilter = fixes.get('headerFilter', False)
fixes.non_ref      = fixes.get('non_ref', False)

if inverse:
	inv = lambda v: not v if isinstance(v, bool) else v
	fixes.clinvarLink  = inv(fixes.clinvarLink)
	fixes.addChr       = inv(fixes.addChr)
	fixes.addAF        = inv(fixes.addAF)
	fixes.tumorpos     = inv(fixes.tumorpos)
	fixes.headerInfo   = inv(fixes.headerInfo)
	fixes.headerContig = inv(fixes.headerContig)
	fixes.headerFormat = inv(fixes.headerFormat)
	fixes.headerFilter = inv(fixes.headerFilter)
	fixes.non_ref      = inv(fixes.non_ref)

if fixes.tumorpos is True:
	fixes.tumorpos = path.splitext(path.basename(infile))[0]
if fixes.tumorpos and isinstance(fixes.tumorpos, str):
	fixes.tumorpos = fixes.tumorpos.split(',')

def getContigsFromFai(fai):
	ret = {}
	with open(fai) as f:
		for line in f:
			parts = line.split('\t')
			ret[parts[0]] = int(parts[1])
	return ret

def getContigsFromDict(dictfile):
	ret = {}
	with open(dictfile) as f:
		for line in f:
			parts = line.split('\t')
			if parts[0] != '@SQ':
				continue
			ret[parts[1][3:]] = int(parts[2][3:])
	return ret

# fix clinvarLink first, because it will fail VCF parser
# just remove it
if fixes.clinvarLink or fixes.addChr or fixes.tumorpos or fixes.non_ref:
	tmpoutfile = outfile + '.tmp'
	openfun = gzip.open if infile.endswith('.gz') else open
	with openfun(infile, 'rt', errors='replace') as fin, \
		openfun(tmpoutfile, 'wt', errors='replace') as fout:
		for line in fin:
			if line.startswith('##'):
				# if no header operated
				fout.write(line)
			elif line.startswith('#') and fixes.tumorpos:
				# determine if we have 2 samples
				parts = line.strip().split('\t')
				if len(parts) != 11:
					fixes.tumorpos = False
				elif parts[10] in fixes.tumorpos:
					fixes.tumorpos = True
				elif parts[9] in fixes.tumorpos:
					fixes.tumorpos = False
				else:
					sample1_match = max(SequenceMatcher(None, parts[9], tumor).ratio()
						for tumor in fixes.tumorpos)
					sample2_match = max(SequenceMatcher(None, parts[10], tumor).ratio()
						for tumor in fixes.tumorpos)
					fixes.tumorpos = sample2_match > sample1_match
				if fixes.tumorpos:
					parts[9], parts[10] = parts[10], parts[9]
				fout.write('\t'.join(parts) + '\n')
			elif line.startswith('#'):
				fout.write(line)
			else:
				parts = line.strip().split('\t')
				if fixes.tumorpos:
					parts[9], parts[10] = parts[10], parts[9]
				if fixes.addChr:
					parts[0] = parts[0] if parts[0].startswith('chr') else 'chr' + parts[0]
				if fixes.non_ref and "<NON_REF>" in parts[4]:
					parts[4] = ','.join(alt for alt in parts[4].split(',') if alt != '<NON_REF>') \
						or '.'
				info = parts[7]
				if fixes.clinvarLink:
					info = ';'.join(inf for inf in info.split(';') if not inf.startswith('<a href'))
				parts[7] = info

				fout.write('\t'.join(parts) + '\n')
	infile = tmpoutfile

pool = {}
if fixes.headerInfo not in (None, False):
	pool['info'] = set()
if fixes.headerContig not in (None, False):
	pool['contig'] = set()
if fixes.headerFormat not in (None, False):
	pool['format'] = set()
if fixes.headerFilter not in (None, False):
	pool['filter'] = set()

if pool or fixes.addAF:
	# scan the problems
	vcf  = VCF(infile)

	for variant in vcf:
		if fixes.headerContig not in (None, False):
			pool['contig'].add(variant.CHROM)
		if fixes.headerInfo not in (None, False):
			for key,_ in variant.INFO:
				pool['info'].add(key)
		if fixes.headerFormat not in (None, False):
			for key in variant.FORMAT:
				pool['format'].add(key)
		if fixes.headerFilter not in (None, False) and variant.FILTER:
			for filt in variant.FILTER.split(';'):
				pool['filter'].add(filt)
	vcf.close() # remove auto-generated headers

	if fixes.headerContig is True:
		fixes.headerContig = {}
	if fixes.headerInfo is True:
		fixes.headerInfo = {}
	if fixes.headerFormat is True:
		fixes.headerFormat = {}
	if fixes.headerFilter is True:
		fixes.headerFilter = {}

	vcf  = VCF(infile)

	# if fixes.addAF is true, try to add it to header
	if fixes.addAF:
		try:
			vcf["AD"]
		except KeyError:
			raise ValueError('Cannot fix addAF, vcf file does not contain "FORMAT/AD" field.')
		try:
			vcf["DP"]
		except KeyError:
			raise ValueError('Cannot fix addAF, Vcf file does not contain "FORMAT/DP" field.')
		try:
			item_type = vcf["AF"]
			if item_type['Description'] == '"Dummy"':
				raise KeyError()
		except KeyError:
			adict = dict(
				ID          = "AF",
				Number      = 1,
				Type        = "Float",
				Description = "The allele frequency of first minor allele."
			)
			vcf.add_format_to_header(adict)

	for info_item in pool.get('info', []):
		try:
			item_type = vcf[info_item]
			# auto added by cyvcf2
			if item_type['Description'] == '"Dummy"':
				raise KeyError()
		except KeyError:
			adict = fixes.headerInfo.get(info_item, {})
			adict = dict(
				ID          = info_item,
				Number      = adict.get('Number', 1),
				Type        = adict.get('Type', "String"),
				Description = adict.get('Description', info_item.upper())
			)
			vcf.add_info_to_header(adict)

	for fmt_item in pool.get('format', []):
		try:
			item_type = vcf[fmt_item]
			if item_type['Description'] == '"Dummy"':
				raise KeyError()
		except KeyError:
			adict = fixes.headerFormat.get(fmt_item, {})
			adict = dict(
				ID          = fmt_item,
				Number      = adict.get('Number', 1),
				Type        = adict.get('Type', "String"),
				Description = adict.get('Description', fmt_item.upper())
			)
			vcf.add_format_to_header(adict)

	for filter_item in pool.get('filter', []):
		try:
			item_type = vcf[filter_item]
		except KeyError:
			desc  = fixes.headerFilter.get(filter_item, filter_item.upper())
			adict = dict(
				ID          = filter_item,
				Description = desc
			)
			vcf.add_filter_to_header(adict)

	contigs2drop = set()
	if pool.get('contig'):
		if path.isfile(reffai):
			refcontigs = getContigsFromFai(reffai)
		elif path.isfile(refdict):
			refcontigs = getContigsFromDict(refdict)
		else:
			refcontigs = {contig_item: 999999999 for contig_item in pool['contig']}

		for contig_item in pool['contig']:
			if contig_item not in vcf.seqnames:
				contig_len = refcontigs.get(contig_item, 999999999)
				vcf.add_to_header("##contig=<ID=%s,length=%s>" % (contig_item, contig_len))

	writer = Writer(outfile, vcf)
	#vcf.close()

	for variant in vcf:
		if variant.CHROM in contigs2drop:
			continue
		if fixes.addAF:
			ADs = variant.format("AD")
			if ADs is None:
				continue
			DPs = variant.format("DP")
			variant.set_format("AF", numpy.array([
				numpy.float(ad[1])/numpy.float(DPs[i][0])
				if DPs[i] and float(DPs[i][0]) > 0 else 0.0
				for i, ad in enumerate(ADs)]))
		writer.write_record(variant)
	writer.close()

else:
	shell.mv(infile, outfile)

with open(errfile, 'r') as f:
	for line in f:
		if line.startswith('[E::vcf_parse_format]'):
			raise ValueError(line)
