{{rimport}}('__init__.r')
options(stringsAsFactors = FALSE)
library(sciClone)
library(vcfR)

mutfile = {{i.vfvcfs | R}}
cnvfile = {{i.cnvcfs | R}}
outdir  = {{o.outdir | quote}}
params  = {{args.params | R}}
exfile  = {{args.exfile | quote}}
mutctrl = {{args.mutctrl | R}}
cnctrl  = {{args.mutctrl | R}}





vfsamcol = {{args.vfsamcol | R}}
cnsamcol = {{args.cnsamcol | R}}
varcount = {{args.varcount | :a if a else 'NULL'}}
cncount  = {{args.cncount | :a if a else 'NULL'}}
if (!is.function(varcount)){
	varcount = function(fmt) as.integer(unlist(strsplit(fmt$AD, ","))[2])
}
if (!is.function(cncount)) {
	cncount = function(fmt) fmt$CN
}
varWithInfo = length(formals(varcount)) > 1
cnWithInfo  = length(formals(cncount)) > 1

if (length(vfvcfs) > 2) {
	stop('No more than 2 samples supported at this moment.')
}

if (length(cnvcfs) > 0 && length(cnvcfs) != length(vfvcfs)) {
	stop('Unequal # copy number data and the input vcf files.')
}

makeSample = function(fmt, sam) {
	# fmt = DP:FDP:SDP:SUBDP:AU:CU:GU:TU
	# sam = 2:0:0:0:0,0:2,3:0,0:0,0	2:0:0:0:2,2:0,0:0,0:0,0
	ns = unlist(strsplit(fmt, ":", fixed = TRUE))
	sm = unlist(strsplit(sam, ":", fixed = TRUE))
	as.list(setNames(sm, ns))
}

makeInfo = function(vari) {
	# SOMATIC;QSS=2;TQSS=1;NT=ref;QSS_NT=2;TQSS_NT=1;SGT=CC->CC;DP=5;MQ=38.23;MQ0=0;ReadPosRankSum=0.00;SNVSB=0.00;SomaticEVS=1.50
	do.call(c, lapply(unlist(strsplit(as.character(vari), ';', fixed = TRUE)), function(x) {
		ret = list()
		if (!grepl('=', x)) {
			ret[[x]] = TRUE
		} else {
			parts = unlist(regmatches(x, regexpr("=", x), invert = TRUE))
			ret[[ parts[1] ]] = parts[2]
		}
		ret
	}))
}

vcf2vaf = function(vcffile) {
	vcf = read.vcfR(vcffile)
	# vcf@gt
	# FORMAT	NORMAL	TUMOR
	# DP:FDP:SDP:SUBDP:AU:CU:GU:TU	2:0:0:0:0,0:2,3:0,0:0,0	2:0:0:0:2,2:0,0:0,0:0,0
	# DP:FDP:SDP:SUBDP:AU:CU:GU:TU	2:0:0:0:0,0:0,0:2,2:0,0	2:0:0:0:2,2:0,0:0,0:0,0
	vaf = data.frame(
		chr      = character(),
		start    = integer(),
		refCount = integer(),
		varCount = integer(),
		VAF      = double()
	)
	samname = colnames(vcf@gt)[vfsamcol + 1]
	for (i in 1:nrow(vcf@gt)) {
		sample = makeSample(vcf@gt[i, "FORMAT"], vcf@gt[i, vfsamcol + 1])
		if (varWithInfo) {
			vari = makeInfo(vcf@fix[i, "INFO"])
			varc = varcount(sample, vari)
		} else {
			varc = varcount(sample)
		}
		if (!is.list(varc)) {
			varc = list(count = varc)
		}
		if (is.null(varc$count)) 
			next
		varc$count = as.numeric(varc$count)
		varc$depth = as.numeric(list.get(varc, 'depth', as.numeric(sample$DP)))
		vaf        = rbind(vaf, list(
			chr      = vcf@fix[i, "CHROM"],
			start    = as.numeric(vcf@fix[i, "POS"]),
			refCount = as.numeric(varc$depth),
			varCount = as.numeric(varc$count),
			VAF      = 100.0*varc$count/varc$depth
		))
	}
	list(vaf = vaf, sample = samname)
}

vcf2cn = function(vcffile) {
	vcf = read.vcfR(vcffile)
	cn  = data.frame(
		chr    = character(),
		start  = integer(),
		stop   = integer(),
		probes = integer(),
		cn     = double()
	)
	#cnfile = file.path(outdir, paste0(vcf@gt[cnsamcol + 1], '.cn'))
	for (i in 1:nrow(vcf@gt)) {
		sample = makeSample(vcf@gt[i, "FORMAT"], vcf@gt[i, cnsamcol + 1])
		cninfo = makeInfo(vcf@fix[i, "INFO"])
		if (cnWithInfo) {
			cno    = cncount(sample, cninfo)
		} else {
			cno = cncount(sample)
		}
		if (!is.list(cno))
			cno = list(cn = cno)
		if (is.null(cno$cn))
			next
		cno$cn     = as.numeric(cno$cn)
		cno$end    = list.get(cno, 'end', cninfo$END)
		cno$probes = list.get(cno, 'probes', cninfo$PROBES)
		cn         = rbind(cn, list(
			chr    = vcf@fix[i, "CHROM"],
			start  = as.numeric(vcf@fix[i, "POS"]),
			stop   = as.numeric(cno$end),
			probes = as.numeric(cno$probes),
			cn     = cno$cn
		))
	}
	cn
}

if (nchar(exfile) > 0 && file.exists(exfile)) {
	exreg = read.table(exfile, header = FALSE, row.names = FALSE, sep = "\t", check.names = FALSE)
} else {
	exreg = NULL
}

vafsams = lapply(vfvcfs, vcf2vaf)
vafs = lapply(vafsams, function(x) x$vaf)
sams = unlist(lapply(vafsams, function(x) x$sample))
rm(vafsams)
if (length(cnvcfs) > 0) {
	cns = lapply(cnvcfs, vcf2cn)
} else {
	cns = NULL
}

params$vafs             = vafs
params$copyNumberCalls  = cns
params$sampleNames      = sams
params$regionsToExclude = exreg

sc = do.call(sciClone, params)
clustfile = file.path(outdir, paste0(paste(sams, collapse = '-'), '.cluster.txt'))
plotfile1 = file.path(outdir, paste0(paste(sams, collapse = '-'), '.plot1d.pdf'))
plotfile2 = file.path(outdir, paste0(paste(sams, collapse = '-'), '.plot2d.pdf'))
writeClusterTable(sc, clustfile)
sc.plot1d(sc, plotfile1)
if (length(sams) > 1) {
	sc.plot2d(sc, plotfile2)
}



