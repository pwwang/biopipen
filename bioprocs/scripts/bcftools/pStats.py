from pyppl import Box
from pathlib import Path
from bioprocs.utils import shell2 as shell

infile        = {{i.infile | quote}}
outdir        = {{o.outdir | quote}}
prefix        = Path(outdir) / {{i.infile | stem | stem | quote}}
bcftools      = {{args.bcftools | quote}}
plot_vcfstats = {{args.plot_vcfstats | quote}}
plot          = {{args.plot | repr}}
params        = {{args.params | repr}}

shell.load_config(
	bcftools = bcftools,
	plot_vcfstats = plot_vcfstats,
)
params._ = infile
params._out = prefix.with_suffix('.stats.txt')
shell.bcftools.stats(**params)

shell.fg.plot_vcfstats(p = outdir, _ = params._out)
