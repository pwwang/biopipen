from pyppl import Box
from bioprocs.utils import shell

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
tool     = {{args.tool | quote}}
gz       = {{args.gz | repr}}
bcftools = {{args.bcftools | quote}}

shell.TOOLS.bcftools = bcftools

def run_bcftools():
	bcftools        = shell.Shell(subcmd = True, equal = ' ').bcftools
	hparams         = Box()
	hparams.h       = True
	hparams._       = infile
	hparams._stdout = outfile
	bcftools.view(**hparams).run()

	qparams = Box()
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
	qparams.f        = "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t.\tGT[\t%GT]\n"
	qparams._        = infile
	qparams._stdout_ = outfile
	if gz:
		qparams.O = 'z'
	bcftools.query(**qparams).run()

tools = dict(
	bcftools = run_bcftools
)

try:
	tools[tool]()
except KeyError:
	raise ValueError('Tool {!r} not supported yet.'.format(tool))
