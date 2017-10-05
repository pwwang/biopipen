from os import path
from tempfile import gettempdir
from pyppl import Proc
from .utils import txt

"""
@name:
	pMTarget2GTargetMat
@description:
	Convert motif target from MSigDB database (i.e. c3.tft.v5.2.entrez.gmt from GSEA to gene-target matrix
	You also have to have a map of motif name to genes (https://github.com/andrewdyates/transcription_factors/blob/master/gsea_msigdb/transfac_id_to_genes_raw.tab)
@input:
	`gmtfile:file`: typically c3.tft.v5.2.entrez.gmt (have to be in entrez format)
	`mapfile:file`: the motif-gene name map file
@output:
	`outfile:file`: the gene-target matrix
@args:
	`species`: The species used to convert gene names, default: human
@requires:
	[python-mygene](https://pypi.python.org/pypi/mygene/3.0.0)
"""
pMTarget2GTargetMat        = Proc ()
pMTarget2GTargetMat.input  = "gmtfile:file, mapfile:file"
pMTarget2GTargetMat.output = "outfile:file:gtmat-{{#}}.txt"
pMTarget2GTargetMat.args   = {'species': 'human'}
pMTarget2GTargetMat.lang   = "python"
pMTarget2GTargetMat.script = """
import mygene
mg = mygene.MyGeneInfo()
# Normalize gene names
# get all genes  in mapfile
genes = {}
origenes = []
with open ("{{mapfile}}") as f:
	for line in f:
		line = line.strip("\\n")
		gene = line.split("\\t")[-1].strip()
		if not gene: continue
		for g in gene.split(';'):
			if g == 'None': continue
			origenes.append(gene)

grets = mg.getgenes (origenes, scopes=['symbol', 'alias'], fields='symbol', species='{{args.species}}')
for gret in grets:
	if not gret.has_key('symbol'): continue
	genes[gret['query']] = gret['symbol']

motif2gene = {}
with open ("{{mapfile}}") as f:
	for line in f:
		line = line.strip("\\n")
		(motif, gene) = line.split("\\t")
		if gene == "None":
			continue
		
		gs = gene.split(";")
		motif2gene[motif] = []
		for g in gs:
			if not genes.has_key(g): continue
			motif2gene[motif].append(genes[g])

# also convert all genes to symbols
allgenes = {}
with open ("{{gmtfile}}") as f:
	for line in f:
		line = line.strip("\\n")
		items = line.split("\\t")
		motif = items.pop(0)
		items.pop(0)
		for item in items:
			allgenes[item] = 1

targets = {}
grets = mg.getgenes (allgenes.keys(), fields='symbol')
for gret in grets:
	if not gret.has_key('symbol'): continue
	targets[gret['query']] = gret['symbol']

final = {}   # tf -> targets
with open ("{{gmtfile}}") as f:
	for line in f:
		line = line.strip("\\n")
		items = line.split("\\t")
		motif = items.pop(0)
		items.pop(0) # url
		if not motif2gene.has_key (motif):
			#sys.stderr.write ("No TF found for motif %s \\n" % (motif))
			continue

		tfs = motif2gene[motif]

		for tf in tfs:
			ori = [] if not final.has_key(tf) else final[tf]
			ori += [targets[g] for g in items if targets.has_key(g)]
			final[tf] = ori


tfs = sorted(final.keys())
allgenes = targets.values()
allgenes = sorted(list(set(allgenes)))
fout = open ("{{outfile}}", "w")
fout.write ("tf\\t" + "\\t".join(tfs) + "\\n")
for g in allgenes:
	fout.write (g)
	for tf in tfs:
		fout.write ("\\t%d" % (1 if g in final[tf] else 0))
	fout.write("\\n")
fout.close()
"""

"""
@name:
	pIntersectGMT
@description:
	Get the intersect gene set from multiple gmt files
	To do intersect for more than 2 files: gmtfile1, gmtfile2, gmtfile3:
	```
	pIntersectGMT.input = {pIntersectGMT.input: channel.create([(gmtfile1, gmtfile2)])}
	
	pIntersectGMT2 = pIntersectGMT.copy()
	pIntersectGMT2.depends = pIntersectGMT
	pIntersectGMT2.input   = {pIntersectGMT2.input.keys()[0]: lambda ch: ch.insert(0, gmtfile3)}
	```
@input:
	`gmtfile1:file`: The 1st gmt file
	`gmtfile2:file`: The 2nd gmt file
@output:
	`outdir:file`: the output gmtfile
@args:
	`geneformat`: The gene names in gene set. Default: "symbol,alias". Available values see mygene docs.
	`gz`: whether the files are with gz format, default: False. If `gz = True`, output file will be also gzipped.
	`species`: The species, used for gene name conversion in mygene, default: human
@requires:
	[python-mygene](https://pypi.python.org/pypi/mygene/3.0.0) 
"""
pIntersectGMT = Proc()
pIntersectGMT.input  = "gmtfile1:file, gmtfile2:file"
pIntersectGMT.output = "outfile:file:intersected.gmt{{args.gz | lambda x: '.gz' if x else ''}}"
pIntersectGMT.args   = {"geneformat": "symbol, alias", "gz":False, "species": "human"}
pIntersectGMT.script = """
#!/usr/bin/env python
import gzip
from mygene import MyGeneInfo
mg = MyGeneInfo()
openfunc = open if not {{args.gz}} else gzip.open
def readGMT (gmtfile):
	ret = {}
	with openfunc (gmtfile) as f:
		for line in f:
			line = line.strip()
			if not line or line.startswith('#'): continue
			items = line.split("\\t")
			ret[items.pop(0)] = [items.pop(1), items]
	return ret
	
gmt1 = readGMT("{{gmtfile1}}")
gmt2 = readGMT("{{gmtfile2}}")
keys = list(set(gmt1) & set(gmt2))

with openfunc ("{{outfile}}", "w") as fout:
	for key in keys:
		comm1 = gmt1[key][0]
		comm2 = gmt2[key][0]
		comm  = comm1+"|"+comm2 if comm1!=comm2 else comm1
		inter = list (set(gmt1[key][1]) & set(gmt2[key][1]))
		genes = mg.querymany (inter, scopes="{{args.geneformat}}", fields="symbol", species="{{args.species}}")
		genes = [gene['symbol'] for gene in genes if gene.has_key('symbol')]
		if not inter: continue
		fout.write ("%s\\t%s\\t%s\\n" % (key, comm, "\\t".join(inter)))
		
"""

"""
@name:
	pUnionGMT
@description:
	Get the union gene set from multiple gmt files
	To do union for more than 2 files: gmtfile1, gmtfile2, gmtfile3:
	```
	pUnionGMT.input = {pIntersectGMT.input: channel.create([(gmtfile1, gmtfile2)])}
	
	pUnionGMT2 = pIntersectGMT.copy()
	pUnionGMT2.depends = pIntersectGMT
	pUnionGMT2.input   = {pUnionGMT.input.keys()[0]: lambda ch: ch.insert(0, gmtfile3)}
	```
@input:
	`gmtfile1:file`: The 1st gmt file
	`gmtfile2:file`: The 2nd gmt file
@output:
	`outdir:file`: the output gmtfile
@args:
	`geneformat`: The gene names in gene set. Default: "symbol,alias". Available values see mygene docs.
	`gz`: whether the files are with gz format, default: False. If `gz = True`, output file will be also gzipped.
	`species`: The species, used for gene name conversion in mygene, default: human
@requires:
	[python-mygene](https://pypi.python.org/pypi/mygene/3.0.0) 
"""
pUnionGMT = Proc()
pUnionGMT.input  = "gmtfile1:file, gmtfile2:file"
pUnionGMT.output = "outfile:file:unioned.gmt{{args.gz | lambda x: '.gz' if x else ''}}"
pUnionGMT.args   = {"geneformat": "symbol, alias", "gz":False, "species": "human"}
pUnionGMT.script = """
#!/usr/bin/env python
import gzip
from mygene import MyGeneInfo
mg = MyGeneInfo()
openfunc = open if not {{args.gz}} else gzip.open
def readGMT (gmtfile):
	ret = {}
	with openfunc (gmtfile) as f:
		for line in f:
			line = line.strip()
			if not line or line.startswith('#'): continue
			items = line.split("\\t")
			ret[items.pop(0)] = [items.pop(1), items]
	return ret
	
gmt1 = readGMT("{{gmtfile1}}")
gmt2 = readGMT("{{gmtfile2}}")
keys = list(set(gmt1) & set(gmt2))

with openfunc ("{{outfile}}", "w") as fout:
	for key in keys:
		comm1 = gmt1[key][0]
		comm2 = gmt2[key][0]
		comm  = comm1+"|"+comm2 if comm1!=comm2 else comm1
		inter = list (set(gmt1[key][1]) | set(gmt2[key][1]))
		genes = mg.querymany (inter, scopes="{{args.geneformat}}", fields="symbol", species="{{args.species}}")
		genes = [gene['symbol'] for gene in genes if gene.has_key('symbol')]
		if not inter: continue
		fout.write ("%s\\t%s\\t%s\\n" % (key, comm, "\\t".join(inter)))
		
"""

"""
@name:
	pExpmat2Gct
@description:
	Convert expression matrix to GCT file.
	Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for file format
@input:
	`expfile:file`: the input expression matrix file. Samples as columns, genes as rows.
@output:
	`outfile:file`: the gct file
"""
pExpmat2Gct        = Proc(desc = 'Convert expression matrix to GCT file.')
pExpmat2Gct.input  = 'expfile:file'
pExpmat2Gct.output = 'outfile:file:{{ in.expfile | fn }}.gct'
pExpmat2Gct.lang   = 'python'
pExpmat2Gct.script = "file:scripts/gsea/pExpmat2Gct.py"

"""
@name:
	pSampleinfo2Cls
@description:
	Convert sample infomation to cls file.
	Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CLS for file format
	NOTE that the order of samples must be the same as in GMT file in further analysis.
@input:
	`sifile:file`: the sample information file.
		- Headers are: [Sample, ]Patient, Group, Batch
		- Rows are samples
@output:
	`outfile:file`: the cls file
"""
pSampleinfo2Cls                       = Proc(desc = 'Convert sample infomation to cls file.')
pSampleinfo2Cls.input                 = 'sifile:file'
pSampleinfo2Cls.output                = 'outfile:file:{{ in.sifile | fn }}.cls'
pSampleinfo2Cls.tplenvs.txtSampleinfo = txt.sampleinfo.python
pSampleinfo2Cls.lang                  = 'python'
pSampleinfo2Cls.script                = "file:scripts/gsea/pSampleinfo2Cls.py"

"""
@name:
	pSSGSEA
@description:
	Single sample GSEA
	Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for GCT file format
	Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GMT for GMT file format
@input:
	`gctfile:file`: the expression file
	`gmtfile:file`: the gmtfile for gene sets
@output:
	`outdir:file`: the output directory
	- `report.txt`: the enrichment report for each Gene set.
	- `RES_<GeneSet>.png`: the running ES plot for <GeneSet>
	- `normP_<GeneSet>.png`: the norminal P value plot for <GeneSet>
@args:
	`weightexp`: Exponential weight employed in calculation of enrichment scores. Default: 0.75
	`nperm`:     Number of permutations. Default: 10000
"""
pSSGSEA = Proc (desc = 'Do single-sample GSEA.')
pSSGSEA.input     = "gctfile:file, gmtfile:file"
pSSGSEA.output    = "outdir:file:{{in.gctfile | fn}}-{{in.gmtfile | fn}}-ssGSEA"
pSSGSEA.args      = {'weightexp': 1, 'nperm': 10000}
pSSGSEA.lang      = 'Rscript'
pSSGSEA.script    = "file:scripts/gsea/pSSGSEA.r"

"""
@name:
	pGSEA
@description:
	GSEA
	Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for GCT file format
	Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GMT for GMT file format
	Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CLS for CLS file format
@input:
	`gctfile:file`: the expression file
	`clsfile:file`: the class file
	`gmtfile:file`: the gmtfile for gene sets
@output:
	`outdir:file`: the output directory
@args:
	`weightexp`: Exponential weight employed in calculation of enrichment scores. Default: 0.75
	`nperm`:     Number of permutations. Default: 10000
"""
pGSEA = Proc (desc = 'Do GSEA.')
pGSEA.input     = "gctfile:file, clsfile:file, gmtfile:file"
pGSEA.output    = "outdir:dir:{{in.gctfile | fn}}-{{in.gmtfile | fn}}-GSEA"
pGSEA.args      = {'weightexp': 1, 'nperm': 1000, 'nthread': 1}
pGSEA.lang      = 'Rscript'
pGSEA.script    = "file:scripts/gsea/pGSEA.r"

"""
@name:
	pEnrichr
@description:
	Use APIs from http://amp.pharm.mssm.edu/Enrichr/help#api&q=1 to analyze a gene list
@input:
	`infile:file`: The gene list, each per line
@output:
	`outdir:dir`:  The output directory, containing the tables and figures.
@args:
	`topn`: Top N pathways used to plot. Default: 10
	`dbs`:  The databases to do enrichment against. Default: KEGG_2016
	  - A full list can be found here: http://amp.pharm.mssm.edu/Enrichr/#stats
	  - Multiple dbs separated by comma (,)
	`norm`: Normalize the gene list use [python-mygene](https://pypi.python.org/pypi/mygene/3.0.0)
	`rmtags`: Remove pathway tags in the plot. Default: True
	  - For example: change "Lysine degradation_Homo sapiens_hsa00310" to "Lysine degradation".
	`plot`: Whether to plot the result. Default: True
	`title`: The title for the plot. Default: "Gene enrichment: {db}"
@requires:
	[python-mygene](https://pypi.python.org/pypi/mygene/3.0.0) if `args.norm` is `True`
"""
pEnrichr              = Proc()
pEnrichr.input        = "infile:file"
pEnrichr.output       = "outdir:dir:{{in.infile | fn}}.enrichr"
pEnrichr.lang         = "python"
pEnrichr.args.topn    = 10
pEnrichr.args.dbs     = "KEGG_2016"
pEnrichr.args.norm    = False
pEnrichr.args.rmtags  = True
pEnrichr.args.plot    = True
pEnrichr.args.title   = "Gene enrichment: {db}"
pEnrichr.args.mgcache = path.join(gettempdir(), 'mygeneinfo.cache')
pEnrichr.beforeCmd    = 'if [[ ! -e "{{args.mgcache}}" ]]; then mkdir -p "{{args.mgcache}}"; fi'
pEnrichr.errhow       = 'retry'
pEnrichr.script       = "file:scripts/gsea/pEnrichr.py"


"""
@name:
	pTargetEnrichr
@description:
	Use APIs from http://amp.pharm.mssm.edu/Enrichr/help#api&q=1 to analyze a gene list
@input:
	`infile:file`: The target genes with regulators
		- Format: 
		- Header is not required, but may specified in first line starting with `#`
		- If only 3 columns are there, the 3rd column is anyway the relation!
		- If only 4 columns are there, 3rd is target status, 4th is relation!
		  ```
		  #Regulator	Target	Regulator status	Target status	Relation
		  has-mir-22	Gene	+	+	+
		  ```
@output:
	`outdir:dir`:  The output directory, containing the tables and figures.
@args:
	`dbs`       : The databases to do enrichment against. Default: KEGG_2016
	  - A full list can be found here: http://amp.pharm.mssm.edu/Enrichr/#stats
	  - Multiple dbs separated by comma (,)
	`rmtags`    : Remove pathway tags in the plot. Default: True
	  - For example: change "Lysine degradation_Homo sapiens_hsa00310" to "Lysine degradation".
	`enrplot`   : Whether to plot the result. Default: True
	`enrn`      : Top N pathways used to plot. Default: 10
	`netplot`   : Whether to plot the network. Default: True
	`netn`      : Top N pathways used to plot the network. Default: 5
		- Must <= `enrn`. If `netn` >= `enrn`, `netn` = `enrn`
	`title`     : The title for the plot. Default: "Gene enrichment: {db}"
@requires:
	[`python-mygene`](https://pypi.python.org/pypi/mygene/3.0.0) 
	[`graphviz`](https://pypi.python.org/pypi/graphviz)
"""
pTargetEnrichr              = Proc(desc = 'Do gene set enrichment analysis for target genes.')
pTargetEnrichr.input        = "infile:file"
pTargetEnrichr.output       = "outdir:dir:{{in.infile | fn}}.tenrichr"
pTargetEnrichr.lang         = "python"
pTargetEnrichr.args.dbs     = "KEGG_2016"
pTargetEnrichr.args.rmtags  = True
pTargetEnrichr.args.enrplot = True
pTargetEnrichr.args.enrn    = 10
pTargetEnrichr.args.netplot = True
pTargetEnrichr.args.netn    = 5
pTargetEnrichr.args.title   = "Gene enrichment: {db}"
pTargetEnrichr.args.mgcache = path.join(gettempdir(), 'mygeneinfo.cache')
pTargetEnrichr.beforeCmd    = 'if [[ ! -e "{{args.mgcache}}" ]]; then mkdir -p "{{args.mgcache}}"; fi'
pTargetEnrichr.errhow       = 'retry'
pTargetEnrichr.script       = "file:scripts/gsea/pTargetEnrichr.py"

