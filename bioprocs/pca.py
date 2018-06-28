from pyppl import Proc
from . import params

"""
@name:
	pPCA
@description:
	Perform PCA analysis
@input:
	`infile:file`: The matrix to do the analysis
	- Note that rows are samples, columns are features, if not, use `args.transpose = True`
@output:
	`outfile:file`: The output coordinate file
	- Columns are PCs, rows are samples
@args:
	`transpose`:  Whether to transpose the input matrix from infile. Default: False
	`rownames`:   The `row.names` argument for `read.table`, default: 1
	`header`:     The `header` argument for `read.table` to read the input file, default: True.
	`screeplot`:  Whether to generate the screeplot or not. Default: True
	`sp_ncp`:     Number of components in screeplot. Default: 0 (auto detect)
	- if total # components (tcp) < 20: use all
	- else if tcp > 20, use 20
	`varplot`:    Whether to generate the variable plot or not. Default: False
	`biplot`:     Whether to generate the variable plot or not. Default: True
@requires:
	[`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html) for plots
"""
pPCA = Proc()
pPCA.input          = "infile:file"
pPCA.output         = "outfile:file:{{infile | fn}}.pca/pcs.txt, outdir:dir:{{infile | fn}}.pca"
pPCA.args.transpose = False
pPCA.args.rownames  = 1
pPCA.args.header    = True
pPCA.args.screeplot = True
pPCA.args.sp_ncp    = 0
pPCA.args.varplot   = False
pPCA.args.biplot    = True
pPCA.lang           = "Rscript"
pPCA.script         = """
library("factoextra")
rownames = {{args.rownames | lambda x: x if str(x).isdigit() else "NULL" }}
data = read.table ("{{infile}}", row.names=rownames, header={{args.header | Rbool }}, check.names=F)
if ({{ args.transpose | Rbool }}) data = t(data)
pc     = prcomp (data)
pcfile = "{{outfile}}"
pcs    = pc$x
tcp    = ncol (pcs)
if ({{args.screeplot | Rbool}}) {
	spfile = file.path ("{{outdir}}", "screeplot.png")
	png (file = spfile)
	ncp = {{args.sp_ncp}}
	if (ncp == 0 && tcp < 20) {
		ncp = tcp
	} else if (ncp == 0 && tcp > 20) {
		ncp = 20
	}
	print (fviz_screeplot (pc, ncp = ncp))
	dev.off()
}
if ({{args.varplot | Rbool}}) {
	vpfile = file.path ("{{outdir}}", "varplot.png")
	png (file = vpfile)
	print (fviz_pca_var(pc, col.var="contrib",
		gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
		repel = TRUE # Avoid text overlapping
	))
	dev.off()
}
if ({{args.biplot | Rbool}}) {
	bpfile = file.path ("{{outdir}}", "biplot.png")
	png (file = bpfile)
	print (fviz_pca_biplot(pc, repel = TRUE))
	dev.off()
}
sdfile = file.path ("{{outdir}}", "stdev.txt")
write.table (pc$sdev, sdfile, col.names=F, row.names=F, quote=F)
write.table (pcs, pcfile, col.names=T, row.names=T, quote=F, sep="\\t")
"""

"""
@name:
	pSelectPCs
@description:
	Select a subset of PCs from pPCA results
@input:
	`indir:file`: The directory generated from pPCA
@output:
	`outfile:file`: The file containing selected PCs
@args:
	`n`: The number of PCs to select. Default: 0.9
	- If it is < 1, used as the % variation explained from stdev.txt
"""
pSelectPCs = Proc ()
pSelectPCs.input  = "infile:file, indir:file"
pSelectPCs.output = "outfile:file:{{indir | fn}}.toppcs.txt"
pSelectPCs.args.n = .9
pSelectPCs.lang   = "python"
pSelectPCs.script = """
import os

n = {{args.n}}
if n < 1:
	stdfile = os.path.join("{{indir}}", "stdev.txt")
	stdevs  = open (stdfile).read().strip().split("\\n")
	stdevs  = map (float, stdevs)
	expstd  = sum (stdevs) * n
	cstd    = 0
	for i, stdev in enumerate(stdevs):
		cstd += stdev
		if cstd > expstd:
			n = i + 1
			break

with open ("{{infile}}") as fin, open ("{{outfile}}", "w") as fout:
	i = 0
	for line in fin:
		line  = line.strip()
		if not line: continue
		parts = line.split("\\t")
		if i == 0:
			parts = parts[:n]
			i += 1
		else:
			parts = parts[:n+1]
		fout.write ("\\t".join(parts) + "\\n")
"""