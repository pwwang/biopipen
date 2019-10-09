from os import path, readlink
from bioprocs.utils import shell2 as shell, gztype

def check(ref):
	if not ref or not path.exists(ref):
		raise Exception('Reference file not exists: %s' % ref)

def checkIndex(refindex):
	if not refindex:
		refindex = [refindex]
	return all([path.exists(ri) for ri in refindex])

def buildIndex(ref, cmd, ref2 = None, cmd2 = None):
	try:
		runcmd(cmd)
		return ref
	except Exception:
		try:
			if not path.exists(ref2):
				shell.ln_s(ref, ref2)
			runcmd(cmd2)
			return ref2
		except Exception:
			return None

def bamIndex(bam, ext = '.bam.bai', samtools = 'samtools', nthread = 1):
	"""
	Index bam files
	If bam file is a link, try to find the index file in its orginal directory or
	its realpath directory
	If nothing found, try to create the index file using samtools
	@params:
		`ext`: The expected extension of index file. Default: `.bam.bai`
			- Some tools requird `XXX.bai` without `.bam`
		`samtools`: The path to samtools. Default: `samtools`
			- If it's None, then an exception will raised instead of creating the index file
		`nthread`: The # threads used to create the index file. Default: `1`
	"""
	if not ext.startswith('.'):
		ext = '.' + ext
	# /path/to/some.bam -> some.bam
	bname = path.basename(bam)
	# /path/to/some.bam -> /path/to/
	dname = path.dirname(bam)
	# some.bam -> some
	fname = path.splitext(bname)[0]
	# some -> some
	# [1]some -> some
	rname = fname.split(']', 1)[1] if fname.startswith('[') else fname

	shell.load_config(samtools = samtools)
	# /path/to/some.bam.bai
	expectedIndex = path.join(dname, fname + ext)
	if path.isfile(expectedIndex):
		return
	# if bam is not a link, there is nowhere else to find index, create it using samtools
	if not path.islink(bam):
		if samtools:
			shell.samtools.index(b = True, _out = expectedIndex, **{'@': nthread})
		else:
			raise ValueError('Index not found: {}'.format(bam))
		return
	# find the index in original directory
	origbam   = readlink(bam)
	origIndex = path.splitext(origbam)[0] + ext
	if path.isfile(origIndex):
		shell.ln_s(origIndex, expectedIndex)
		return
	# find the index in realpath directory
	realbam   = path.realpath(bam)
	realIndex = path.splitext(realbam)[0] + ext
	if path.isfile(realIndex):
		shell.ln_s(realIndex, expectedIndex)
		return
	# if all failed, create it
	if samtools:
		shell.samtools.index(b = True, _out = expectedIndex, **{'@': nthread})
	else:
		raise ValueError('Index not found: {}'.format(bam))

def tabixIndex(filename, type, tabix = 'tabix'):
	# /path/to/some.vcf -> some.vcf
	# /path/to/some.vcf.gz -> some.vcf
	bname = path.basename(filename[:-3]) if filename.endswith('.gz') else path.basename(filename)
	# /path/to/some.bam -> /path/to/
	dname = path.dirname(filename)
	# some.vcf -> some
	# some.vcf.gz -> some
	fname = path.splitext(bname[:-3] if bname.endswith('.gz') else bname)[0]
	# some -> some
	# [1]some -> some
	rname = fname.split(']', 1)[1] if fname.startswith('[') else fname

	expectedIndex = path.join(dname, fname + '.%s.gz.tbi' % type)
	if path.isfile(expectedIndex):
		return filename if filename.endswith('.gz') else filename + '.gz'
	
	shell.load_config(tabix = tabix)
	# type could bed3, bed6..
	ptype = 'bed' if type.startswith('bed') else type
	gt = gztype(filename)
	if gt == 'bgzip':
		if path.islink(filename):
			linkfile = readline(filename)
			if path.isfile(linkfile + '.tbi'):
				shell.ln_s(linkfile + '.tbi', expectedIndex)
				return filename
			readfile = path.realpath(filename)
			if path.isfile(realfile + '.tbi'):
				shell.ln_s(realfile + '.tbi', expectedIndex)
				return realfile
		shell.fg.tabix(p = ptype, _ = filename)
		return filename
	if gt == 'gzip':
		bgzfile = path.join(dname, bname + '.bgz.' + type)
		shell.gunzip(filename, c = True, _out = bgzfile)
		shell.bgzip(bgzfile)
		shell.fg.tabix(p = ptype, _ = bgzfile)
		return bgzfile + '.gz'
	shell.bgzip(filename, c = True, _out = filename + '.gz')
	shell.fg.tabix(p = ptype, _ = filename + '.gz')
	return filename + '.gz'

def vcfIndex(vcf, tabix = 'tabix'):
	return tabixIndex(vcf, 'vcf', tabix)

def bedIndex(bed, tabix = 'tabix'):
	return tabixIndex(bed, 'bed', tabix)
