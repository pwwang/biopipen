from hashlib import sha256
from os import path
from diot import Diot
from bioprocs.utils import shell2 as shell

seed         = {{i.seed | repr}}
fq1          = {{o.fq1 | quote}}
fq2          = {{o.fq2 | quote}}
tool         = {{args.tool | quote}}
wgsim        = {{args.wgsim | quote}}
dwgsim       = {{args.dwgsim | quote}}
art_illumina = {{args.art_illumina | quote}}
len1         = {{args.len1 | repr}}
len2         = {{args.len2 | repr}}
num          = {{args.num | repr}}
gz           = {{args.gz | repr}}
params       = {{args.params | repr}}
ref          = {{args.ref | quote}}

shell.load_config(wgsim = wgsim, dwgsim = dwgsim, art_illumina = art_illumina)

def get_seed(s):
	if s is None:
		return -1
	if isinstance(s, int):
		return s
	return int(sha256(str(s).encode()).hexdigest()[:7], 16)

def run_wgsim():
	global fq1, fq2
	if gz:
		fq1, fq2 = fq1[:-3], fq2[:-3]
	params.N = num
	params['1'] = len1
	params['2'] = len2
	params.S = get_seed(seed)
	params._ = [ref, fq1, fq2]
	shell.wgsim(**params).fg
	if gz:
		shell.gzip(fq1)
		shell.gzip(fq2)

def run_dwgsim():
	global fq1, fq2
	if gz:
		fq1, fq2 = fq1[:-3], fq2[:-3]
	prefix = path.commonprefix([fq1, fq2]).rstrip('._[]')

	params.N = num
	params['1'] = len1
	params['2'] = len2
	params.S = 2
	params.z = get_seed(seed)
	params._ = [ref, prefix]
	shell.dwgsim(**params).fg

	shell.mv(prefix + '.bwa.read1.fastq', fq1)
	shell.mv(prefix + '.bwa.read2.fastq', fq2)

	if gz:
		shell.gzip(fq1)
		shell.gzip(fq2)

def run_art_illumina():
	global fq1, fq2
	if gz:
		fq1, fq2 = fq1[:-3], fq2[:-3]
	prefix = path.commonprefix([fq1, fq2]).rstrip('._[]')

	params.c = num
	params.i = ref
	params.l = len1
	params.o = prefix
	params.p = True
	params.rndSeed = get_seed(seed)

	shell.art_illumina(**params).fg

	shell.mv(prefix + '1.fq', fq1)
	shell.mv(prefix + '2.fq', fq2)

	if gz:
		shell.gzip(fq1)
		shell.gzip(fq2)

tools = dict(
	wgsim = run_wgsim,
	dwgsim = run_dwgsim,
	art_illumina = run_art_illumina
)

try:
	tools[tool]()
except KeyError:
	raise KeyError('Tool %r not supported.' % tool)
