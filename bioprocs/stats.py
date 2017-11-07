from pyppl import Proc
from . import params
from .utils import helpers

"""
@name:
	pMetaPval
@description:
	Combine p-values in the files from input directory
@input:
	`indir:dir`: The directory containing the input files
@output:
	`outfile:file`: The output file containing the meta-pvalues
@args:
	`args.pattern`: The pattern used to filter the input files. Default: '*'
	`args.header`: Whether the input files contains a header. Default: True
		- Could be a list to specify it for each file.
		- The order should be concordant with the file names
	`args.pcol`: Which column is the p-value. Default: -1 (last column)
	`args.poutonly`: Only output pvalues. Default: False (output all possible information)
	`args.outheader`: Whether output the header. Default: True
	`args.method`: The method used to calculate the meta-pvalue. Default: sumlog (Fisher's method)
		- Other available methods: logitp, sumz, votep, sump, meanp and wilkinsonp
		- See: https://www.rdocumentation.org/packages/metap/versions/0.8
@requires:
	[`r-matep`](https://www.rdocumentation.org/packages/metap/)
"""
pMetaPval                   = Proc(desc = "Combine p-values.")
pMetaPval.input             = "indir:dir"
pMetaPval.output            = "outfile:file:{{in.indir | fn}}.metapval.txt"
pMetaPval.args.pattern      = "*"
pMetaPval.args.header       = True
pMetaPval.args.pcol         = -1
pMetaPval.args.poutonly     = False
pMetaPval.args.outheader    = True
pMetaPval.args.method       = 'sumlog' # fisher's method
pMetaPval.tplenvs.cbindfill = helpers.cbindfill.r
pMetaPval.lang              = params.Rscript.value
pMetaPval.script            = "file:scripts/stats/pMetaPval.r"

"""
@name:
	pMetaPval1
@description:
	Combine p-values in a single file by rows.
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output file containing the meta-pvalues
@args:
	`args.header`: Whether the input files contains a header. Default: True
	`args.pcol`: Which column is the p-value. Default: -1 (last column)
	`args.poutonly`: Only output pvalues. Default: False (output all possible information)
	`args.outheader`: Whether output the header. Default: True
	`args.method`: The method used to calculate the meta-pvalue. Default: sumlog (Fisher's method)
		- Other available methods: logitp, sumz, votep, sump, meanp and wilkinsonp
		- See: https://www.rdocumentation.org/packages/metap/versions/0.8
@requires:
	[`r-matep`](https://www.rdocumentation.org/packages/metap/)
"""
pMetaPval1 = Proc(desc = "Combine p-values in a single file by rows.")
pMetaPval1.input             = "infile:file"
pMetaPval1.output            = "outfile:file:{{in.infile | fn}}.metapval.txt"
pMetaPval1.args.header       = True
pMetaPval1.args.pcol         = -1
pMetaPval1.args.poutonly     = False
pMetaPval1.args.outheader    = True
pMetaPval1.args.method       = 'sumlog' # fisher's method
pMetaPval1.lang              = params.Rscript.value
pMetaPval1.script            = "file:scripts/stats/pMetaPval1.r"