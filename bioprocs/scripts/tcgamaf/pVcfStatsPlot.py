######### from csvStats to matrix files
from os import path, makedirs
from math import ceil
from glob import glob
import sys

{{runcmd}}

chroms  = []
chparts = "{{args.chroms}}".split(',')
chparts = map (lambda x: x.strip(), chparts)
for chpart in chparts:
	if '-' not in chpart:
		chroms.append (chpart)
		continue
	ch1, ch2 = chpart.split('-')
	ch1 = int(ch1)
	ch2 = int(ch2)
	chroms += map(str, list (range(ch1, ch2+1)))
chroms = filter(None, chroms)

matxdir = "{{o.outdir}}/mats"
if not path.isdir(matxdir): makedirs(matxdir)
plotdir = "{{o.outdir}}/plots"
if not path.isdir(plotdir): makedirs(plotdir)
plotsrc = "{{o.outdir}}/plot.r"

flags  = {}
data   = {}
def flagHit (name, sampe):
	if not flags.has_key(name):
		flags[name] = {}
		return False
	if not flags[name].has_key(sample): return False
	return flags[name][sample]

def dataHit (name, sample):
	if not data.has_key(name):
		data[name] = {}
		return False
	if not data[name].has_key(sample): return False
	return True

def handleFreq (sample, line, name, starts, ends, skips, datacol = 1):
	if not line.startswith (starts) and not flagHit(name, sample):
		return
	
	if line.startswith (starts):
		print "- Start: %s" % name
		if not flagHit (name, sample):
			flags[name][sample] = True
		if not dataHit (name, sample):
			data[name][sample]  = {}
		return
	
	if skips:
		if not isinstance(skips, list): skips = skips.split(",")
		for skip in skips:
			if line.startswith (skip): return
		
	if line.startswith (ends):
		flags[name][sample] = False
		return
	
	parts  = line.replace("%", "").split(",")
	key    = parts[0].strip()
	dat    = parts[datacol].strip()
	if name in ['Changes_by_chromosome', 'Change_rate_by_chromosome'] and chroms and key not in chroms:
		return
	
	data[name][sample][key] = dat

	
def handleScatters (sample, line, name, starts, ends, keys, vals, datacol = 1):
	if not line.startswith (starts) and not flagHit(name, sample):
		return
		
	if line.startswith (starts):
		print "- Start: %s" % name
		if not flagHit (name, sample):
			flags[name][sample] = True
		if not dataHit (name, sample):
			data[name][sample]  = {}
		return
		
	if line.startswith (ends):
		flags[name][sample] = False
		return
		
	parts  = line.replace("%", "").split(",")
	parts  = map(lambda x: x.strip(), parts)
	if line.startswith (keys):
		data[name][sample]['keys'] = parts[datacol:]
	if line.startswith (vals):
		data[name][sample]['vals'] = parts[datacol:]

def handleChangeTable (sample, line, name, starts, ends, headstr):
	if not line.startswith (starts) and not flagHit(name, sample):
		return
		
	if line.startswith (starts):
		print "- Start: %s" % name
		if not flagHit (name, sample):
			flags[name][sample] = True
		if not dataHit (name, sample):
			data[name][sample]  = {}
		return
		
	if line.startswith (ends):
		flags[name][sample] = False
		del data[name][sample]["head"]
		return
		
	if line.startswith (headstr):
		parts  = line.split(",")
		parts  = map(lambda x: x.strip(), parts)
		parts.pop(0)
		data[name][sample]["head"] = parts
		return
	
	parts  = line.replace("%", "").split(",")
	parts  = map(lambda x: x.strip(), parts)
	k      = parts.pop(0)
	if not data[name][sample].has_key(k): data[name][sample][k] = {}
	for i, h in enumerate(data[name][sample]["head"]):
		data[name][sample][k][h] = parts[i]
		
def handleChromChange (sample, line, name, starts, datacol = 2):
	if not line.startswith (starts) and not flagHit(name, sample):
		return
		
	if line.startswith (starts):
		print "- Start: %s" % name
		if not flagHit (name, sample):
			flags[name][sample] = True
		if not data.has_key(name):
			data[name] = {}
		return
	
	if "Position" in line or "Count" in line:
		parts  = line.split(",")
		parts  = map(lambda x: x.strip(), parts)
		chrom  = parts[0]
		if chroms and chrom not in chroms:
			return
		if not dataHit (name, chrom): data[name][chrom] = {}
		if not data[name][chrom].has_key(sample):
			data[name][chrom][sample] = {}
		data[name][chrom][sample][parts[1]] = parts[datacol:]
		
def saveFreqHist (name, sampledata):
	outfile  = path.join(matxdir, name + '.mat.txt')
	# get keys/colnames
	colnames = []
	for sample, sdata in sampledata.iteritems(): colnames += sdata.keys()
	colnames = list(set(colnames))
	colnames.sort(key=lambda x: float(x) if x.isdigit() else x)
	
	fout    = open (outfile, 'w')
	fout.write ("\t" + "\t".join (colnames) + "\n")
	for sample, sdata in sampledata.iteritems():
		fout.write (sample + "\t" + "\t".join (['0' if not sdata.has_key(cname) else sdata[cname] for cname in colnames]) + "\n")
	fout.close()
	
def saveScatters (name, sampledata):
	outfile  = path.join(matxdir, name + '.mat.txt')
	sums     = {}
	colnames = []
	for sample, sdata in sampledata.iteritems():
		if not sdata.has_key('keys'): continue
		colnames += sdata['keys']
		for i, key in enumerate(sdata['keys']):
			val = sdata['vals'][i]
			if not sums.has_key(key): sums[key] = [0, 0]
			sums [key][0] += float(val)
			sums [key][1] += 1
	means    = {key:val[0]/val[1] for key, val in sums.iteritems()}
	
	colnames = list(set(sums.keys()))
	colnames.sort(key=lambda x: float(x) if x.isdigit() else x)
	if not colnames: return
	fout     = open (outfile, 'w')
	fout.write ("\t" + "\t".join(colnames) + "\n")
	
	for sample, sdata in sampledata.iteritems():
		outs  = [sample]
		outs += [means[key] if key not in sdata['keys'] else sdata['vals'][sdata['keys'].index(key)] for key in colnames]
		fout.write ("\t".join(outs) + "\n")
	fout.close()
	
def saveChangeTable (name, sampledata):
	outdir   = path.join (matxdir, name)
	makedirs (outdir)
	names    = []
	for sample, sdata in sampledata.iteritems():
		names += sdata.keys()
		if names:
			names += sdata[names[0]].keys()
	names = list(set(names))
	names.sort(key=lambda x: float(x) if x.isdigit() else x)
	
	for sample, sdata in sampledata.iteritems():
		outfile = path.join (outdir, sample + ".mat.txt")
		fout    = open (outfile, 'w')
		fout.write ("\t" + "\t".join(names) + "\n")
		for rname in names:
			outs    = [rname]
			if sdata.has_key (rname):
				outs += ['0' if not sdata[rname].has_key(cname) else sdata[rname][cname] for cname in names]
			else:
				outs += ['0'] * len (names)
			fout.write ("\t".join(outs) + "\n")
		fout.close()
		
def saveChromChange (name, sampledata):
	outdir   = path.join (matxdir, name)
	makedirs (outdir)
	for chrom, sdata in sampledata.iteritems():
		outfile = path.join (outdir, chrom + ".mat.txt")
		colnames= []
		fout    = open (outfile, 'w')
		for sample, dat in sdata.iteritems():
			if not colnames:
				colnames = dat['Position']
				fout.write ("\t" + "\t".join(dat['Position']) + "\n")
			elif colnames != dat['Position']:
				sys.stderr.write ("Warning: chomosome %s bins are different" % chrom)
				sys.exit (1)
			fout.write (sample + "\t" + "\t".join(dat['Count']) + "\n")
		fout.close()

for i, stfile in enumerate(glob("{{i.indir}}/*.stats*.csv")):
	sample                             = path.basename(stfile).split('.')[0]
	print "Handling %s ..." % stfile
	f = open (stfile)
	for line in f:
		line = line.strip()
		if not line: continue
		handleFreq    (sample, line, "Change_rate", "# Summary table", "# Change rate by chromosome", "N,G,D,S,Co,W", 1)
		handleFreq    (sample, line, "Changes_by_chromosome", "# Change rate by chromosome", "# ", "Chromosome", 2)
		handleFreq    (sample, line, "Change_rate_by_chromosome", "# Change rate by chromosome", "# ", "Chromosome", 3)
		handleFreq    (sample, line, "Effect_count_by_impact", "# Effects by impact", "# Effects by functional class", "Type", 1)
		handleFreq    (sample, line, "Effect_percentage_by_impact", "# Effects by impact", "# Effects by functional class", "Type", 2)
		handleFreq    (sample, line, "Effect_count_by_functional_class", "# Effects by functional class", "# Count by effects", "Type,Missense_Silent_ratio", 1)
		handleFreq    (sample, line, "Effect_percentage_by_functional_class", "# Effects by functional class", "# Count by effects", "Type,Missense_Silent_ratio", 2)
		handleFreq    (sample, line, "Missense_Silent_ratio", "# Effects by functional class", "# Count by effects", "Type,MISSENSE,NONSENSE,SILENT", 1)
		handleFreq    (sample, line, "Count_by_effects", "# Count by effects", "# Count by genomic region", "Type", 1)
		handleFreq    (sample, line, "Percentage_by_effects", "# Count by effects", "# Count by genomic region", "Type", 2)
		handleFreq    (sample, line, "Count_by_genomic_region", "# Count by genomic region", "# Quality", "Type", 1)
		handleFreq    (sample, line, "Percentage_by_genomic_region", "# Count by genomic region", "# Quality", "Type", 2)
		handleScatters    (sample, line, "Indel_length", "# InDel lengths", "# Base changes", "Values", "Count", 1)
		handleChangeTable (sample, line, "Base_change", "# Base changes", "# Ts/Tv summary", "base")
		handleFreq    (sample, line, "TsTv_count", "# Ts/Tv summary", "# Ts/Tv : All variants", "Ts_Tv_ratio", 1)
		handleFreq    (sample, line, "TsTv_ratio", "# Ts/Tv summary", "# Ts/Tv : All variants", "Transitions,Transversions", 1)
		handleFreq    (sample, line, "HomHet_count", "# Hom/Het table", "# Codon change table", "Sample_names", 1)
		handleChangeTable (sample, line, "Condon_change", "# Codon change table", "# Amino acid change table", "codons")
		handleChangeTable (sample, line, "Amino_acid_change", "# Amino acid change table", "# Chromosome change table", "aa")
		handleChromChange (sample, line, "Chromosome_change", "# Chromosome change table")
	
	flags["Chromosome_change"][sample] = False
	f.close()

for name, sampledata in data.iteritems():
	if name in ['Change_rate', 'Changes_by_chromosome', 'Change_rate_by_chromosome', 'Effect_count_by_impact',
				'Effect_percentage_by_impact', 'Effect_count_by_functional_class', 'Effect_percentage_by_functional_class',
				'Missense_Silent_ratio', 'Count_by_effects', 'Percentage_by_effects', 'Count_by_genomic_region',
				'Percentage_by_genomic_region', 'TsTv_count', 'TsTv_ratio', 'HomHet_count']:
		saveFreqHist (name, sampledata)
		continue
		
	if name in ['Indel_length']:
		saveScatters (name, sampledata)
		continue
	
	if name in ['Base_change', 'Condon_change', 'Amino_acid_change']:
		saveChangeTable (name, sampledata)
		continue
	
	if name in ['Chromosome_change']:
		saveChromChange (name, sampledata)
		continue

plotSource = """
library(ggplot2)
library(corrplot)
matxdir = "%s"
plotdir = "%s"

{{plotHist}}
{{plotBoxplot}}

plotFreq = function (fn, xlab="", ylab="Frequency") {
	infile  = file.path (matxdir, paste(fn, ".mat.txt", sep=""))

	tryCatch ({
		outfile = file.path (plotdir, paste(fn, ".png", sep=""))
		data    = read.table (infile, header=T, row.names=1, sep="\\t", check.names=F, stringsAsFactors=FALSE)
		data    = as.matrix(data[,1])

		ggs     = {{args.histplotggs | Rlist}}
		ggs     = c(list(ggtitle(gsub("_", " ", fn)), xlab(xlab), ylab(ylab)), ggs)
		plotHist(data, outfile, devpars={{args.devpars | Rlist}}, ggs=ggs)
	}, error = function (e) {
		write(paste("Failed to plot", fn, ':', e), stderr())
	})
}

plotBoxPlot = function (fn, pie=TRUE, xlab="", ylab="Frequency", shortx = "", shorty = "", marbtm=8) {
	infile  = file.path (matxdir, paste(fn, ".mat.txt", sep=""))
	tryCatch ({
		outfile = file.path (plotdir, paste(fn, ".png", sep=""))
		data    = read.table (infile, header=T, row.names=1, sep="\\t", check.names=F, stringsAsFactors=FALSE)
		data    = as.matrix(data)

		ggs     = {{args.boxplotggs | Rlist}}
		ggs     = c(list(ggtitle(gsub("_", " ", fn)), xlab(xlab), ylab(ylab)), ggs)
		plotBoxplot(data, outfile, devpars={{args.devpars | Rlist}}, ggs=ggs)

		if (pie) {
			samples = rownames(data)
			labels  = colnames(data)
			sampdir = file.path (plotdir, fn)
			dir.create(sampdir, showWarnings = F, recursive = T)
			for (i in 1:nrow(data)) {
				sample = samples[i]
				samplefile = file.path (sampdir, paste(sample, ".png", sep=""))
				png (samplefile, res=300, width=2000, height=2000)
				pie (data[i, ], labels=labels, main=paste(gsub("_", " ", fn), ": ", sample, sep=""))
				dev.off()
			}
		}
	}, error = function (e) {
		write(paste("Failed to plot", fn, ':', e), stderr())
	})
}

plotScatter = function (fn) {
	infile  = file.path (matxdir, paste(fn, ".mat.txt", sep=""))
	tryCatch ({
		if (file.exists(infile)) {
			outfile = file.path (plotdir, paste(fn, ".png", sep=""))
			png (outfile, res=300, width=2000, height=2000)
			data    = read.table (infile, header=T, row.names=1, sep="\\t", check.names=F, stringsAsFactors=FALSE)
			means   = colMeans(data)
			sds     = apply(data, 2, sd)
			maxnt   = 20
			maxt    = ncol(data)
			x       = 1:maxt
			nt      = floor (maxt/ceiling(maxt/maxnt))
			ticks   = seq (1, maxt, by = nt)
			xlabs   = colnames(data)[ticks]
			xlabs   = as.numeric (xlabs) / 1000000
			xlabs   = paste (xlabs, "M", sep="")
			p       = qplot(x, means) +
					  geom_errorbar(aes(x=x, ymin=means-sds, ymax=means+sds), width=2, colour="#999999") +
					  theme_get () +
					  geom_line(col="#333333") +
					  ggtitle (paste(gsub("_", " ", fn), unlist(strsplit(infile, ".", fixed=T))[1], sep=": ")) +
					  labs(x="Chromosome region", y="Changes") +
					  scale_x_continuous (breaks = ticks, labels = xlabs) +
					  theme(axis.text.x = element_text(angle=60, vjust=0.5, hjust=0.5))
			print (p)
			dev.off()
		}
	}, error = function (e) {
		write(paste("Failed to plot", fn, ':', e), stderr())
	})
}

plotChangeTable = function (fn) {
	indir = file.path (matxdir, fn)
	setwd (indir)
	mats  = list()
	sums  = NULL
	files = list.files()
	tryCatch ({
		for (infile in files) {
			mats[[infile]] = as.matrix(read.table (infile, header=T, row.names=1, check.names=F, sep="\\t"))
			if (is.null(sums)) {
				sums = mats[[infile]]
			} else {
				sums = sums + mats[[infile]]
			}
		}
		means = sums / length(files)
		sdsum = NULL
		for (infile in files) {
			if (is.null(sdsum)) {
				sdsum = (mats[[infile]] - means)^2
			} else {
				sdsum = sdsum + (mats[[infile]] - means)^2
			}
		}
		sds   = sqrt (sdsum/(length(files) - 1))
		outfile = file.path (plotdir, paste(fn, ".png", sep=""))
		png (outfile, res=300, width=2000, height=2000)
		diag(means) = 0
		diag(sds)   = 0
		corrplot (means, sizes=sds, is.corr=F, tl.srt = 60, tl.cex=.6, tl.col="black", cl.cex=.6, main=gsub("_", " ", fn), mar=c(1,1,3,1))
		dev.off()
	}, error = function (e) {
		write(paste("Failed to plot", fn, ':', e), stderr())
	})
}

plotChromChange = function (fn) {
	indir   = file.path (matxdir, fn)
	setwd (indir)
	for (infile in list.files()) {
		tryCatch ({
			outfile = file.path (plotdir, paste(fn, "_", unlist(strsplit(infile, ".", fixed=T))[1], ".png", sep=""))
			png (outfile, res=300, width=2000, height=2000)
			data    = as.matrix(read.table (infile, header=T, row.names=1, sep="\\t", check.names=F, stringsAsFactors=FALSE))
			means   = colMeans(data)
			sds     = apply(data, 2, sd)
			maxnt   = 20
			maxt    = ncol(data)
			x       = 1:maxt
			nt      = floor (maxt/ceiling(maxt/maxnt))
			ticks   = seq (1, maxt, by = nt)
			xlabs   = as.numeric(colnames(data)[ticks])
			if (max(xlabs) < 1000000) {
				xlabs = paste (xlabs/1000, "K", sep="")
			} else {
				xlabs   = paste (xlabs/1000000, "M", sep="")
			}
			p       = qplot(x, means) +
					  geom_errorbar(aes(x=x, ymin=means-sds, ymax=means+sds), width=2, colour="#999999") +
					  theme_get () +
					  geom_line(col="#333333") +
					  ggtitle (paste(gsub("_", " ", fn), unlist(strsplit(infile, ".", fixed=T))[1], sep=": ")) +
					  labs(x="Chromosome region", y="Changes") +
					  scale_x_continuous (breaks = ticks, labels = xlabs) +
					  theme(axis.text.x = element_text(angle=60, vjust=0.5, hjust=0.5))
			print (p)
			dev.off()
		}, error = function (e) {
			write(paste("Failed to plot", fn, ':', e), stderr())
		})
	}
}

plotFreq ("Change_rate", xlab="Change rate (1 change/N bps)")
plotBoxPlot ("Changes_by_chromosome", ylab="Changes", xlab="Chromosome")
plotBoxPlot ("Change_rate_by_chromosome", pie=F, ylab="Change rate (1 change/N bps)", xlab="Chromosome")
plotBoxPlot ("Effect_count_by_impact")
plotBoxPlot ("Effect_percentage_by_impact", pie=F, ylab="Percentage")
plotBoxPlot ("Effect_count_by_functional_class")
plotBoxPlot ("Effect_percentage_by_functional_class", pie=F, ylab="Percentage")
plotFreq ("Missense_Silent_ratio")
plotBoxPlot ("Count_by_effects")
plotBoxPlot ("Percentage_by_effects", pie=F, ylab="Percentage")
plotBoxPlot ("Count_by_genomic_region")
plotBoxPlot ("Percentage_by_genomic_region", pie=F, ylab="Percentage")
plotScatter ("Indel_length")
plotChangeTable ("Base_change")
plotChangeTable ("Condon_change")
plotBoxPlot ("TsTv_count")
plotFreq ("TsTv_ratio")
plotBoxPlot ("HomHet_count")
plotChangeTable ("Amino_acid_change")
plotChromChange ("Chromosome_change")
""" % (matxdir, plotdir)

with open(plotsrc, 'w') as f: f.write(plotSource)
runcmd('{{args.Rscript}} "%s"' % plotsrc)
