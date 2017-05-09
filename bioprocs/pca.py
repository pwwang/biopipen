from pyppl import proc

"""
@name:
	pPCA
@description:
	Perform PCA analysis
@input:
	`infile:file`: The matrix to do the analysis
	- Note that rows are features, columns are samples, if not, use `args.transpose = True`
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
pPCA = proc()
pPCA.input   = "infile:file"
pPCA.output  = "outdir:dir:{{infile.fn}}.pca"
pPCA.args    = {
	'transpose': False,
	'rownames': 1,
	'header': True,
	'screeplot': True,
	'sp_ncp': 0,
	'varplot': False,
	'biplot': True
}
pPCA.lang    = "Rscript"
pPCA.script  = """
library("factoextra")
rownames = {{proc.args.rownames | (lambda x: x if str(x).isdigit() else "NULL" )(_)}}
data = read.table ("{{infile}}", row.names=rownames, header={{proc.args.header | str(_).upper() }}, check.names=F)
if ({{ proc.args.transpose | str(_).upper() }}) data = t(data)
pc     = prcomp (data)
pcfile = file.path("{{outdir}}", "pcs.txt")
pcs    = pc$rotation
tcp    = ncol (pcs)
if ({{proc.args.screeplot | str(_).upper()}}) {
	spfile = file.path ("{{outdir}}", "screeplot.png")
	png (file = spfile)
	ncp = {{proc.args.sp_ncp}}
	if (ncp == 0 && tcp < 20) {
		ncp = tcp
	} else if (ncp == 0 && tcp > 20) {
		ncp = 20
	}
	print (fviz_screeplot (pc, ncp = ncp))
	dev.off()
}
if ({{proc.args.varplot | str(_).upper()}}) {
	vpfile = file.path ("{{outdir}}", "varplot.png")
	png (file = vpfile)
	print (fviz_pca_var(pc, col.var="contrib",
		gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
		repel = TRUE # Avoid text overlapping
	))
	dev.off()
}
if ({{proc.args.biplot | str(_).upper()}}) {
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
pSelectPCs = proc ()
pSelectPCs.input  = "indir:file"
pSelectPCs.output = "outfile:file:{{indir.fn}}.pcs.txt"
pSelectPCs.args   = {"n": 0.9}
pSelectPCs.lang   = "python"
pSelectPCs.script = """
import os

n = {{proc.args.n}}
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

with open (os.path.join("{{indir}}", "pcs.txt")) as fin, open ("{{outfile}}", "w") as fout:
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