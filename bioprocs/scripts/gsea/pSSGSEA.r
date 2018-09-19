# Adopted from GSEA R package (GSEA.Gct2Frame2), read the GCT file
readGCT <- function(filename = "NULL") {
      content <- readLines(filename)
      content <- content[-1]
      content <- content[-1]
      col.names <- noquote(unlist(strsplit(content[1], "\t")))
      col.names <- col.names[c(-1, -2)]
      num.cols <- length(col.names)
      content <- content[-1]
      num.lines <- length(content)


      row.nam <- vector(length=num.lines, mode="character")
      #row.des <- vector(length=num.lines, mode="character")
      m <- matrix(0, nrow=num.lines, ncol=num.cols)

      for (i in 1:num.lines) {
         line.list <- noquote(unlist(strsplit(content[i], "\t")))
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
		line = unlist(strsplit(temp[[i]], "\t"))
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
		qs = p.adjust (pvec)
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
	write.table (outmat, outfile, sep= "\t", quote=F, col.names=T, row.names=T)
	#write.table (es, outfile, append=T, sep= "\t", quote=F, col.names=T, row.names=T)
}
seed = {{args.seed}}
if (seed > -1) set.seed(seed)
dir.create("{{o.outdir}}", showWarnings = F, recursive = T)
exp = readGCT("{{i.gctfile}}")
gs  = readGMT("{{i.gmtfile}}")
ess = ESWithPerm(exp = exp, gs = gs)
esq = AdjustP(ess)
ESPlot(es = esq, gs = gs, outprefix = "{{o.outdir}}/RES.")
nes = NPvalPlot(es = esq, outprefix = "{{o.outdir}}/normP.")
ExportResult (es = esq, nes=nes, outfile = "{{o.outdir}}/report.txt")
