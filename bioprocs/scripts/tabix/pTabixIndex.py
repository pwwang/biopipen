import binascii, sys
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

def gzip_type(fn):
	with open(fn, 'rb') as f:
		flag = binascii.hexlify(f.read(4))
	if flag == b'1f8b0804':
		return 'bgzip'
	if flag == b'1f8b0808':
		return 'gzip'
	return 'flat'

def gunzip(fn, outfn = None):
	args = {'f': True, '1': True}
	if not outfn:
		cmd = 'gunzip {!r} {}'.format(fn, args)
	else:
		cmd = 'gunzip {!r} {} -c > {!r}'.format(fn, args, outfn)
	runcmd(cmd)

def bgzip(fn, outfn = None):
	if not outfn:
		cmd = 'bgzip {!r}'.format(fn)
	else:
		cmd = 'bgzip {!r} -c > {!r}'.format(fn, outfn)
	runcmd(cmd)

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}

gz     = {{args.gz | repr}}
params = {{args.params | repr}}
tabix  = {{args.tabix | quote}}

# make sure outfile without .gz
if gz: outfile = outfile[:-3]

gztype = gzip_type(infile)

if gztype == 'gzip':
	gunzip(infile, outfile)
	if gz:
		bgzip(outfile)
elif gztype == 'bgzip':
	if gz:
		sys.symlink(infile, outfile + '.gz')
	else:
		gunzip(infile, outfile)
else:
	if gz:
		bgzip(infile, outfile + '.gz')
	else:
		sys.symlink(infile, outfile)

cmd = '{} {} {!r}'.format(tabix, cmdargs({'p': 'vcf'}), outfile + '.gz' if gz else outfile)
runcmd(cmd)
