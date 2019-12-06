{{rimport}}('__init__.r')
library(QuantumClone)
library(vcfR)
options(stringsAsFactors = FALSE)

vfvcfs   = {{i.vfvcfs | R}}
outdir   = {{o.outdir | quote}}
params   = {{args.params | R}}
vfsamcol = {{args.vfsamcol | R}}
varcount = {{args.varcount | :a if a else 'NULL'}}
nthread  = {{args.nthread | R}}
if (!is.function(varcount)){
	varcount = function(fmt) as.integer(unlist(strsplit(fmt$AD, ","))[2])
}
varWithInfo = length(formals(varcount)) > 1

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
		SampleName = character(),
		Chr        = character(),
		Start      = integer(),
		Depth      = integer(),
		Alt        = integer(),
		Genotype   = character()
	)
	samname = colnames(vcf@gt)[vfsamcol + 1]
	for (i in 1:nrow(vcf@gt)) {
		sample = makeSample(vcf@gt[i, "FORMAT"], vcf@gt[i, vfsamcol + 1])
		tryCatch({
			if (varWithInfo) {
				vari = makeInfo(vcf@fix[i, "INFO"])
				varc = varcount(sample, vari)
			} else {
				varc = varcount(sample)
			}
		}, error = function(e) {
			logger('Unable to get AD of variant ', vcf@fix[i, "CHROM"], ':', vcf@fix[i, "POS"], level = 'WARNING')
			varc = NULL
		})
		if (!is.list(varc)) {
			varc = list(count = varc)
		}
		if (is.null(varc$count)) 
			next
		varc$count = as.numeric(varc$count)
		varc$depth = as.numeric(list.get(varc, 'depth', as.numeric(sample$DP)))
		if (sample$GT == '.') {
			logger('Unknown genotype, use AB instead for', vcf@fix[i, "CHROM"], ':', vcf@fix[i, "POS"], level = 'WARNING')
			varGt = 'AB'
		} else {
			gt1 = as.integer(substr(sample$GT, 1, 1))
			gt2 = as.integer(substr(sample$GT, nchar(sample$GT), nchar(sample$GT)))
			if (gt1 + gt2 == 0) {
				varGt = 'AA'
			} else if (gt1 + gt2 == 2) {
				varGt = 'BB'
			} else {varGt = 'AB'}
		}
		vaf        = rbind(vaf, list(
			SampleName = samname,
			Chr        = vcf@fix[i, "CHROM"],
			Start      = as.numeric(vcf@fix[i, "POS"]),
			Depth      = as.numeric(varc$depth),
			Alt        = as.numeric(varc$count),
			Genotype   = varGt
		))
	}
	rownames(vaf) = NULL
	vaf[complete.cases(vaf),,drop = FALSE]
}

# prepare input data
indata = list()
for (vfvcf in vfvcfs) {
	indata = c(indata, list(vcf2vaf(vfvcf)))
}

output = One_step_clustering(
	indata, 
	contamination    = rep(0, length(vfvcfs)),
	output_directory = outdir
)

