from pyppl import Box
from bioprocs.utils import shell
from bioprocs.utils.parallel import Parallel
from bioprocs.utils.reference import vcfIndex

invcfs    = {{i.infiles | repr}}
outfile   = {{o.outfile | quote}}
nthread   = {{args.nthread | int}}
joboutdir = {{job.outdir | quote}}
vcftools  = {{args.vcftools | quote}}
gatk      = {{args.gatk | quote}}
tabix     = {{args.tabix | quote}}
ref       = {{args.ref | quote}}
params    = {{args.params | repr}}
tool      = {{args.tool | quote}}

shell.TOOLS.vcftools = vcftools
shell.TOOLS.gatk     = gatk

para = Parallel(nthread, raiseExc = True)
invcfs = para.run(vcfIndex, [
	(vcf, tabix) for vcf in invcfs
])

def run_vcftools():
	params.d       = params.get('d', True)
	params.t       = params.get('t', True)
	params._       = invcfs
	params._stdout = outfile
	shell.Shell(equal = ' ').vcftools(**params).run()

def run_gatk():
	params.T                   = 'CombineVariants'
	params.o                   = outfile
	params.R                   = ref
	params.nt                  = nthread
	params.variant             = invcfs
	params._stdout             = outfile
	params.genotypemergeoption = params.get('genotypemergeoption', 'UNIQUIFY')
	shell.Shell(equal = ' ', duplistkeys = True).gatk(**params).run()

tools = dict(
	vcftools = run_vcftools,
	gatk     = run_gatk
)

try:
	tools[tool]()
except KeyError:
	raise KeyError('Tool {!r} not supported.'.format(tool))
except:
	raise
finally:
	pass