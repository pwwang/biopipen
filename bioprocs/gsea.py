from pyppl import proc

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
pMTarget2GTargetMat = proc ()
pMTarget2GTargetMat.input     = "gmtfile:file, mapfile:file"
pMTarget2GTargetMat.output    = "outfile:file:gtmat.txt"
pMTarget2GTargetMat.args      = {'species': 'human'}
pMTarget2GTargetMat.defaultSh = "python"
pMTarget2GTargetMat.script    = """
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

grets = mg.getgenes (origenes, scopes=['symbol', 'alias'], fields='symbol', species='{{proc.args.species}}')
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
pSSGSEA = proc ()
pSSGSEA.input     = "gctfile:file, gmtfile:file"
pSSGSEA.output    = "outdir:file:{{gctfile.fn}}_results"
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
EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = {{proc.args.weightexp}}, correl.vector = NULL) {  
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
		if (maxRES < 0.3) maxRES <- 0.3
		if (minRES > -0.3) minRES <- -0.3
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
		sub.string <- paste("ES =", signif(esrets[[1]]$ES, digits = 3), ", NES =", signif(nes, digits=3), ", Nom. p-val=", signif(esrets[[1]]$p, digits = 3), ", FDR=", signif(esrets[[1]]$q, digits = 3), sep="", collapse="")

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

ESWithPerm = function (exp, gs, nperm = {{proc.args.nperm}}) {
  exp      = exp[order(-exp[,1]),,drop=F]
  corrVec  = c (exp[,1])
  
  geneList = rownames(exp)
  ret = list()
  for (gset in names(gs)) {
    ess = vector (length = nperm, mode = "numeric")
    ret[[ gset ]] = vector (length = nperm, mode = "list")
    ret[[ gset ]][[1]] = EnrichmentScore (geneList, gs[[gset]], correl.vector = corrVec)
    ess[1] = ret[[ gset ]][[1]]$ES
    for (i in 1:nperm-1) {
      geneList1 = sample(geneList)
      ret[[ gset ]][[i+1]] = EnrichmentScore (geneList1, gs[[gset]], correl.vector = corrVec)
      ess[i+1] = ret[[ gset ]][[i+1]]$ES
    }
    ess.sorted = sort (ess, index.return = T)
    for (i in 1:nperm) {
      ret[[ gset ]][[i]]$p = ess.sorted$ix[i] / nperm
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
		qs = p.adjust (pvec, "{{proc.args.padjust}}")
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
	}
	outmat = data.frame (cbind(ES = ES, NES = nes, Pval = pval, Qval = fdr, PeatAt = maxI))
	rownames (outmat) = names(es)
	write.table (outmat, outfile, sep= "\\t", quote=F, col.names=T, row.names=T)
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


