{{rimport}}('plot.r', 'sampleinfo.r')

library(methods)
library(simpleaffy)

indir     = {{i.indir | R}}
saminfo   = {{i.sifile | R}}
pattern   = {{args.pattern | R}}
norm      = {{args.norm | R}}
cdffile   = {{args.cdffile | R}}
outfile   = {{o.outfile | R}}
celdir    = file.path({{job.outdir | quote}}, 'CELs')
prefix    = {{i.indir, args.pattern | *dirpat2name | R}}
transfm   = {{args.transfm | lambda x: x or 'NULL'}}
fn2sample = Vectorize({{args.fn2sample}})

dir.create(celdir, showWarnings = TRUE, recursive = TRUE)
setwd(celdir)
files = Sys.glob(file.path(indir, pattern))
for (celfile in files) {
	tryCatch({
		file.link(celfile, file.path(celdir, basename(celfile)))
	}, error = function(e){
		file.copy(celfile, file.path(celdir, basename(celfile)), overwrite = TRUE)
	}, warning = function(w){
		file.copy(celfile, file.path(celdir, basename(celfile)), overwrite = TRUE)
	})
}

files = basename(files)
samples = fn2sample(files)

gfile = "covdesc.txt"
if (saminfo == '') {
	# generate a covdesc file
	group = matrix(1, ncol = 1, nrow = length(files))
	colnames(group) = c("treatment")
	rownames(group) = files
	write.table(group, gfile, quote = F, sep = "\t")
} else {
	sampleinfo = SampleInfo(saminfo)$dataframe()
	group = matrix(NA, ncol = 1, nrow = length(files))
	colnames(group) = c("treatment")
	rownames(group) = files
	for (i in 1:length(files)) {
		group[files[i], "treatment"] = sampleinfo$Group[which(rownames(sampleinfo) == samples[i])]
	}
	write.table(group, gfile, quote = F, sep = "\t")
}

affydata = read.affy(covdesc = gfile)
if (cdffile != '') {
  library(makecdfenv)
  cdfname = cleancdfname(whatcdf(files[1]))
  assign(cdfname, make.cdf.env(basename(cdffile), cdf.path = dirname(cdffile)))
  affydata@cdfName = cdfname
}
exprs    = call.exprs(affydata, algorithm = norm, do.log = TRUE)
exprsout = exprs@assayData$exprs
colnames(exprsout) = fn2sample(colnames(exprsout))

if (!is.null(transfm)) exprsout = transfm(exprsout)
write.table(round(exprsout, 3), outfile, quote=F, sep="\t")
