from diot import Diot
from bioprocs.utils import shell2 as shell

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
tool     = {{args.tool | quote}}
gz       = {{args.gz | repr}}
bcftools = {{args.bcftools | quote}}

shell.load_config(bcftools = bcftools)

def run_bcftools():
	# bcftools        = shell.Shell(subcmd = True, equal = ' ').bcftools
	hparams         = Diot()
	hparams.h       = True
	hparams._       = infile
	# hparams._stdout = outfile
	shell.bcftools.view(**hparams).r > outfile

	qparams = Diot()
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
	qparams.f        = "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t.\tGT[\t%GT]\n"
	qparams._        = infile
	# qparams._stdout_ = outfile
	if gz:
		qparams.O = 'z'
	shell.bcftools.query(**qparams).r > outfile

tools = dict(
	bcftools = run_bcftools
)

try:
	tools[tool]()
except KeyError:
	raise ValueError('Tool {!r} not supported yet.'.format(tool))
