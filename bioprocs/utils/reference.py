from os import path, symlink, readlink
from bioprocs.utils import runcmd, cmdargs
from bioprocs.utils import shell

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
				symlink(ref, ref2)
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
	
	samtools = Shell({'samtools': samtools}, subcommand = True).samtools if samtools else None
	# /path/to/some.bam.bai 
	expectedIndex = path.join(dname, rname + ext)
	if path.isfile(expectedIndex):
		return
	# if bam is not a link, there is nowhere else to find index, create it using samtools
	if not path.islink(bam):
		if samtools:
			samtools.index(b = True, _stdout = expectedIndex, **{'@': nthread})
		else:
			raise ValueError('Index not found: {}'.format(bam))
		return
	# find the index in original directory
	origbam   = readlink(bam)
	origIndex = path.splitext(origbam)[0] + ext
	if path.isfile(origIndex):
		symlink(origIndex, expectedIndex)
		return
	# find the index in realpath directory
	realbam   = path.realpath(bam)
	realIndex = path.splitext(realbam)[0] + ext
	if path.isfile(realIndex):
		symlink(realIndex, expectedIndex)
		return
	# if all failed, create it
	if samtools:
		samtools.index(b = True, _stdout = expectedIndex, **{'@': nthread})
	else:
		raise ValueError('Index not found: {}'.format(bam))



	
	

