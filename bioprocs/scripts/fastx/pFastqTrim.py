from pyppl import Box
from bioprocs.utils import mem2, shell2 as shell
from pathlib import Path

fq1         = {{i.fq1 | quote}}
fq2         = {{i.fq2 | quote}}
outfq1      = {{o.outfq1 | quote}}
outfq2      = {{o.outfq2 | quote}}
tool        = {{args.tool | quote}}
cutadapt    = {{args.cutadapt | quote}}
skewer      = {{args.skewer |quote}}
trimmomatic = {{args.trimmomatic | quote}}
params      = {{args.params | repr}}
nthread     = {{args.nthread | repr}}
gz          = {{args.gz | repr}}
mem         = {{args.mem | quote }}
minlen      = {{args.minlen |repr}}
minq        = {{args.minq |repr}}
cut5        = {{args.cut5 |repr}}
cut3        = {{args.cut3 |repr}}
adapter1    = {{args.adapter1 |quote}}
adapter2    = {{args.adapter2 |quote}}
jobout      = Path({{job.outdir |quote}})

shell.load_config(
	cutadapt    = cutadapt,
	skewer      = skewer,
	trimmomatic = trimmomatic,
)

def _seqrev (seq):
	d = {
		'A':'T', 'T':'A', 'G':'C', 'C':'G',
		'a':'t', 't':'a', 'g':'c', 'c':'g'
	}
	return ''.join([d[s] for s in seq])

def run_trimmomatic():
	global mem, minlen
	mem    = mem2(mem, 'java')
	minlen = minlen * 2
	adfile = jobout / "adapters.fa"
	adfile.write_text("\n".join([
		">PE1", _seqrev(adapter1),
		">PE1_rc", adapter1,
		">PE2", _seqrev(adapter2),
		">PE2_rc", adapter2,
	]))
	params.nthread = nthread
	params._prefix = '-'
	shell.fg.trimmomatic(
		mem, "PE", params, fq1, fq2, outfq1, "/dev/null", outfq2, "/dev/null",
		"ILLUMINACLIP:%s:2:30:10" % adfile,
		"LEADING:%s" % cut5,
		"TRAILING:%s" % cut3,
		"SLIDINGWINDOW:4:%s" % minq,
		"MINLEN:%s" % minlen
	)

def run_cutadapt():
	params.a = adapter1
	params.A = adapter2
	params.u = [cut5, -cut3]
	params.U = [cut5, -cut3]
	params.m = minlen
	params.q = "{0}.{0}".format(minq)
	params.o = outfq1
	params.O = outfq2
	params._ = [fq1, fq2]
	shell.fg.cutadapt(**params)

def run_skewer():
	params.m = 'pe'
	params.t = nthread
	params.x = adapter1
	params.y = adapter2
	params.Q = minq
	params.l = minlen
	params.z = gz
	params.o = jobout / 'tmp'
	params._ = [fq1, fq2]
	shell.fg.skewer(**params)
	out1 = jobout / 'tmp-trimmed-pair1.fastq'
	out2 = jobout / 'tmp-trimmed-pair2.fastq'
	if gz:
		out1 = str(out1) + '.gz'
		out2 = str(out2) + '.gz'
	shell.mv(out1, outfq1)
	shell.mv(out2, outfq2)

TOOLS = dict(
	trimmomatic = run_trimmomatic,
	cutadapt    = run_cutadapt,
	skewer      = run_skewer
)

try:
	TOOLS[tool]()
except KeyError:
	raise ValueError('{} not supported yet.'.format(tool))
