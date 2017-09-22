from pyppl import Proc

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
	pSSGSEA
@description:
	Single sample GSEA
	Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for file format
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
	`padjust`:   P value adjustment method, default 'bonferroni'. Can be "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
	`nperm`:     Number of permutations. Default: 10000
@requires:
	[python-mygene](https://pypi.python.org/pypi/mygene/3.0.0)
"""
pSSGSEA = Proc ()
pSSGSEA.input     = "gctfile:file, gmtfile:file"
pSSGSEA.output    = "outdir:file:{{gctfile|fn}}_results"
pSSGSEA.args      = {'weightexp': 0.75, 'padjust': 'bonferroni', 'nperm': 10000}
pSSGSEA.script    = """
#!/usr/bin/env Rscript
# Adopted from GSEA R package (GSEA.Gct2Frame2), read the GCT file
readGCT <- function(filename = "NULL") { 
      content <- readLines(filename)
      content <- content[-1]
      content <- content[-1]
      col.names <- noquote(unlist(strsplit(content[1], "\\t")))
      col.names <- col.names[c(-1, -2)]
      num.cols <- length(col.names)
      content <- content[-1]
      num.lines <- length(content)


      row.nam <- vector(length=num.lines, mode="character")
      #row.des <- vector(length=num.lines, mode="character")
      m <- matrix(0, nrow=num.lines, ncol=num.cols)

      for (i in 1:num.lines) {
         line.list <- noquote(unlist(strsplit(content[i], "\\t")))
         row.nam[i] <- noquote(line.list[1])
         #row.des[i] <- noquote(line.list[2])
         line.list <- line.list[c(-1, -2)]
         for (j in 1:length(line.list)) {
            m[i, j] <- as.numeric(line.list[j])
         }
      }
      ds <- data.frame(m)
      names(ds) <- col.names
      row.names(ds) <- row.nam
      return(ds)
}

readGMT <- function (filename = "NULL") {
	temp = readLines (filename)
	max.Ng <- length(temp)
	gs = list()
	for (i in 1:max.Ng) {
		line = unlist(strsplit(temp[[i]], "\\t"))
		name = make.names(line[1])
		line = line [-1]
		line = line [-1] # desc
		gs[[name]] = line
		
	} 
	
	return (gs)
}

# Calculate ES score, adopted from GSEA R package (GSEA.EnrichmentScore)
EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = {{args.weightexp}}, correl.vector = NULL) {  
	tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
	no.tag.indicator <- 1 - tag.indicator 
	N <- length(gene.list) 
	Nh <- length(gene.set) 
	Nm <-  N - Nh 
	if (weighted.score.type == 0) {
		correl.vector <- rep(1, N)
	}
	alpha <- weighted.score.type
	correl.vector <- abs(correl.vector**alpha)
	sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
	norm.tag    <- 1.0/sum.correl.tag
	norm.no.tag <- 1.0/Nm
	RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
	max.ES <- max(RES)
	min.ES <- min(RES)
	if (is.na(max.ES) || is.na(min.ES)) {
		return(list(ES = 0, arg.ES = 1, RES = RES, indicator = tag.indicator))   
	}
	if (max.ES > - min.ES) {
	#      ES <- max.ES
		ES <- signif(max.ES, digits = 5)
		arg.ES <- which.max(RES)
	} else {
	#      ES <- min.ES
		ES <- signif(min.ES, digits=5)
		arg.ES <- which.min(RES)
	}
	return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}

ESPlot = function (es, gs, outprefix) {
	for (gset in names(es)) {
		png (filename = paste (outprefix, make.names(gset), ".png", sep="", collapse = ""), type="cairo")
		esret = es [[gset]][[1]]
		N = length(esret$indicator)
		sub.string <- paste("Number of genes: ", N, " (in list), ", length(gs[[gset]]), " (in gene set)", sep = "", collapse="")

		main.string <- paste("Gene Set: ", gset)
		minRES = min (esret$RES)
		maxRES = max (esret$RES)
		if (is.na(maxRES) || maxRES < 0.3) maxRES <- 0.3
		if (is.na(minRES) || minRES > -0.3) minRES <- -0.3
		delta <- (maxRES - minRES)*.5
		minplot <- minRES - delta
		maxplot <- maxRES
		col <- ifelse(esret$ES > 0, 2, 4)
		plot(1:N, esret$RES, main = main.string, sub = sub.string, xlab = "Gene List Index", ylab = "Running Enrichment Score (RES)", xlim=c(1, N), ylim=c(minplot, maxplot), type = "l", lwd = 2, cex = 1, col = col)

		rect (0, minplot, N, minplot + 0.5*delta, col=rgb(255, 250, 230, maxColorValue=255), border=NA, lwd=0)
		lines(c(1, N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
		lines(c(esret$arg.ES, esret$arg.ES), c(minplot, maxplot), lwd = 1, lty = 3, cex = 1, col = col) # max enrichment vertical line
		for (j in 1:N) {
			if (esret$indicator[j] == 1) {
				lines(c(j, j), c(minplot, minplot + 0.5*delta), lwd = 1, lty = 1, cex = 1, col = 1)  # enrichment tags
			}
		}


		adjx <- ifelse(esret$ES > 0, 0, 1)

		leg.txt <- paste("Peak at ", esret$arg.ES, sep="", collapse="")
		text(x=esret$arg.ES, y=minplot + 1.8*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)
		dev.off()
	}
}

NPvalPlot = function (es, outprefix) {
	Ns = length(es[[1]])
	NESret = vector (length = 0, mode = "numeric")
	for (gset in names(es)) {
		png (filename = paste (outprefix, make.names(gset), ".png", sep="", collapse = ""), type="cairo")
		esrets = es [[gset]]
		nes = esrets[[1]]$ES / mean (esrets[[1]]$RES[esrets[[1]]$RES > 0])
		NESret = c (NESret, nes)
		if (esrets[[1]]$ES < 0) {
			nes = - esrets[[1]]$ES / mean (esrets[[1]]$RES[esrets[[1]]$RES < 0])
		}
		#sub.string <- paste("ES =", signif(esrets[[1]]$ES, digits = 3), ", NES =", signif(nes, digits=3), ", Nom. p-val=", signif(esrets[[1]]$p, digits = 3), ", FDR=", signif(esrets[[1]]$q, digits = 3), sep="", collapse="")
		sub.string <- paste("ES =", signif(esrets[[1]]$ES, digits = 3), ", NES =", signif(nes, digits=3), ", Nom. p-val=", signif(esrets[[1]]$p, digits = 3), sep="", collapse="")

		phi = vector (length = Ns, mode = "numeric")
		for (i in 1:Ns) {
			phi[i] = esrets[[i]]$ES
		}

		temp <- density(phi, adjust=.5)
		x.plot.range <- range(temp$x)
		y.plot.range <- c(-0.125*max(temp$y), 1.5*max(temp$y))
		plot(temp$x, temp$y, type = "l", sub = sub.string, xlim = x.plot.range, ylim = y.plot.range, lwd = 2, col = 2, main = "Gene Set Null Distribution", xlab = "ES", ylab="P(ES)")
		x.loc <- which.min(abs(temp$x - esrets[[1]]$ES))
		lines(c(esrets[[1]]$ES, esrets[[1]]$ES), c(0, temp$y[x.loc]), lwd = 2, lty = 1, cex = 1, col = 1)
		lines(x.plot.range, c(0, 0), lwd = 1, lty = 1, cex = 1, col = 1)

		leg.txt <- c("Gene Set Null Density", "Observed Gene Set ES value")
		c.vec <- c(2, 1)
		lty.vec <- c(1, 1)
		lwd.vec <- c(2, 2)
		legend(x=-0.2, y=y.plot.range[2], bty="n", bg = "white", legend=leg.txt, lty = lty.vec, lwd = lwd.vec, col = c.vec, cex = 1.0)
		dev.off()
	}
	return (NESret)
}

ESWithPerm = function (exp, gs, nperm = {{args.nperm}}) {
  exp      = exp[order(-exp[,1]),,drop=F]
  corrVec  = c (exp[,1])
  
  geneList = rownames(exp)
  ret = list()
  for (gset in names(gs)) {
    ess = vector (length = nperm, mode = "numeric")
    ret[[ gset ]] = vector (length = nperm, mode = "list")
    ret[[ gset ]][[1]] = EnrichmentScore (geneList, gs[[gset]], correl.vector = corrVec)
    ret[[ gset ]][[1]]$hits = geneList[as.logical(ret[[ gset ]][[1]]$indicator)]
    ess[1] = ret[[ gset ]][[1]]$ES
    for (i in 1:(nperm-1)) {
      geneList1 = sample(geneList)
      ret[[ gset ]][[i+1]] = EnrichmentScore (geneList1, gs[[gset]], correl.vector = corrVec)
      #ret[[ gset ]][[i+1]]$hits = geneList1[as.logical(ret[[ gset ]][[i+1]]$tag.indicator)]
      ess[i+1] = ret[[ gset ]][[i+1]]$ES
    }
    #ess.sorted = sort (ess, index.return = T)
    for (i in 1:nperm) {
      if (ess[i] > 0) {
        ret[[ gset ]][[i]]$p = sum (ess >= ess[i]) / nperm
      } else {
        ret[[ gset ]][[i]]$p = sum (ess <= ess[i]) / nperm
      }
    }
  }
  
  return (ret)
}

AdjustP = function (esp) {
	ret = esp
	for (gset in names(ret)) {
		pvec = vector (length = length(ret[[ gset ]]), mode="numeric")
		for (j in 1:length(ret[[ gset ]])) {
			pvec[j] = ret[[ gset ]][[ j ]]$p
		}
		qs = p.adjust (pvec, "{{args.padjust}}")
		for (j in 1:length(ret[[ gset ]])) {
			ret[[ gset ]][[ j ]]$q = qs[j]
		}
	}
	return (ret)
}

ExportResult = function (es, nes, outfile) {
	ES   = vector(length = 0, mode = "numeric")
	pval = vector(length = 0, mode = "numeric")
	fdr  = vector(length = 0, mode = "numeric")
	maxI = integer(length = 0)
	
	for (gset in names(es)) {
		esret = es[[gset]][[1]]
		ES    = c (ES, esret$ES)
		pval  = c (pval, esret$p)
		fdr   = c (fdr, esret$q)
		maxI  = c (maxI, esret$arg.ES)
		Hits  = paste(esret$hits, collapse=",")
	}
	outmat = data.frame (cbind(ES = ES, NES = nes, Pval = pval, Qval = fdr, PeatAt = maxI, Hits=Hits))
	rownames (outmat) = names(es)
	write.table (outmat, outfile, sep= "\\t", quote=F, col.names=T, row.names=T)
	#write.table (es, outfile, append=T, sep= "\\t", quote=F, col.names=T, row.names=T)
}
dir.create("{{outdir}}", showWarnings = F, recursive = T)
exp = readGCT("{{gctfile}}")
gs  = readGMT("{{gmtfile}}")
ess = ESWithPerm(exp = exp, gs = gs)
esq = AdjustP(ess)
ESPlot(es = esq, gs = gs, outprefix = "{{outdir}}/RES.")
nes = NPvalPlot(es = esq, outprefix = "{{outdir}}/normP.")
ExportResult (es = esq, nes=nes, outfile = "{{outdir}}/report.txt")

"""

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
	  - if `topn` < 1: use it as a p-value, otherwise use it as a number cutoff
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
pEnrichr = Proc()
pEnrichr.input  = "infile:file"
pEnrichr.output = "outdir:dir:{{in.infile | fn}}.enrichr"
pEnrichr.lang   = "python"
pEnrichr.args   = {"topn": 10, "dbs": "KEGG_2016", "norm": False, "rmtags": True, "plot": True, "title": "Gene enrichment: {db}"}
pEnrichr.errhow = 'retry'
pEnrichr.script = """
import json
import requests
import math
genes = [line.split()[0] for line in open("{{in.infile}}") if line.strip()]
if {{args.norm}}:
	from mygene import MyGeneInfo
	mg = MyGeneInfo()
	genes = mg.getgenes (genes, scopes=['symbol', 'alias'], fields='symbol', species='human')
	genes = [gene['symbol'] for gene in genes if gene.has_key('symbol')]
	
if {{args.plot}}:
	import matplotlib
	matplotlib.use('Agg')
	from matplotlib import pyplot as plt
	from matplotlib import gridspec
	from matplotlib import patches

## upload
ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
genes_str   = "\\n".join(genes)
description = '{{in.infile | fn}}'
payload = {
    'list': (None, genes_str),
    'description': (None, description)
}

response = requests.post(ENRICHR_URL, files=payload)
if not response.ok:
    raise Exception('Error analyzing gene list')

data = json.loads(response.text)

## do enrichment
dbs = "{{args.dbs}}".split(',')
dbs = map (lambda s: s.strip(), dbs)

ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
query_string = '?userListId=%s&backgroundType=%s'

head = ["#Rank", "Term name", "P-value", "Z-score", "Combined score", "Overlapping genes", "Adjusted p-value", "Old p-value", "Old adjusted p-value"]
topn = {{args.topn}}
for db in dbs:
	user_list_id = data['userListId']
	gene_set_library = db
	response = requests.get(
		ENRICHR_URL + query_string % (user_list_id, gene_set_library)
	)
	if not response.ok:
		raise Exception('Error fetching enrichment results against %s' % db)
	
	data    = json.loads(response.text)
	data    = data[db]
	d2plot  = []
	outfile = "{{out.outdir}}/%s.txt" % db
	fout    = open (outfile, "w")
	fout.write ("\\t".join(head) + "\\n")
	for i, ret in enumerate(data):
		fout.write ("\\t".join(['|'.join(r) if x == 5 else str(r) for x,r in enumerate(ret)]) + "\\n")
		if topn < 1 and ret[2] >= topn: continue
		if topn >= 1 and i > topn - 1: continue
		if {{args.rmtags}} and "_Homo sapiens_hsa" in ret[1]: ret[1] = ret[1][:-22]
		d2plot.append (ret)
	fout.close()
	
	if {{args.plot}}:
		#d2plot   = sorted (d2plot, cmp=lambda x,y: 0 if x[2] == y[2] else (-1 if x[2] < y[2] else 1))
		plotfile = "{{out.outdir}}/%s.png" % db
		gs = gridspec.GridSpec(1, 2, width_ratios=[3, 7]) 
		rownames = [r[1] for r in d2plot]
		rnidx    = range (len (rownames))
		ax1 = plt.subplot(gs[0])
		plt.title ("{{args.title}}".replace("{db}", db), fontweight='bold')
		
		ax1.xaxis.grid(alpha=.6, ls = '--', zorder = -99)
		plt.subplots_adjust(wspace=.01, left=0.5)
		ax1.barh(rnidx, [len(r[5]) for r in d2plot], color='blue', alpha=.6)
		plt.yticks (rnidx, rownames)
		ax1.yaxis.set_ticks_position('none')
		ax1.tick_params(axis='x', colors='blue')
		ax1.spines['top'].set_visible(False)
		ax1.spines['left'].set_visible(False)
		ax1.spines['right'].set_visible(False)
		ax1.spines['bottom'].set_linewidth(1)
		ax1.invert_xaxis()
		ax1.invert_yaxis()
		xticks = ax1.xaxis.get_major_ticks()
		xticks[0].label1.set_visible(False)
		
		ax2 = plt.subplot(gs[1])
		ax2.xaxis.grid(alpha=.6, ls = '--', zorder = -99)
		ax2.barh(rnidx, [-math.log(r[2], 10) for r in d2plot], color='red', alpha = .6)
		for i, r in enumerate(d2plot):
			t  = str("%.2E" % r[2])
			tx = 0.1
			ty = i + 0.1
			ax2.text(tx, ty, t, fontsize=8)
		ax2.tick_params(axis='x', colors='red')
		ax2.spines['top'].set_visible(False)
		ax2.spines['left'].set_visible(False)
		ax2.spines['right'].set_visible(False)
		ax2.spines['bottom'].set_linewidth(1)
		ax2.invert_yaxis()
		plt.yticks([])
		ng_patch = patches.Patch(color='blue', alpha=.6, label='# overlapped genes')
		pv_patch = patches.Patch(color='red', alpha = .6, label='-log(p-value)')
		plt.figlegend(handles=[ng_patch, pv_patch], labels=['# overlapped genes', '-log(p-value)'], loc="lower center", ncol=2, edgecolor="none")
		plt.savefig(plotfile, dpi=300)
"""

