library(methods)
{{rimport}}('__init__.r')

indir  = {{i.indir | R}}
pdir   = {{i.pdir | R}}
outdir = {{o.outdir | R}}
plink  = {{args.plink | R}}

bedfile = Sys.glob(file.path(indir, '*.bed'))
input   = tools::file_path_sans_ext(bedfile)
output  = file.path(outdir, basename(input))

if (dir.exists(pdir)) {
	fails = Sys.glob(file.path(pdir, '*.fail'))
} else {
	fails = pdir
}

isSampleFail = function(fail) {
	endsWith(fail, '.samplecr.fail') ||
	endsWith(fail, '.het.fail') ||
	endsWith(fail, '.ibd.fail')
}

isSNPFail = function(fail) {
	endsWith(fail, '.snpcr.fail') ||
	endsWith(fail, '.hwe.fail')
}

for (i in c('.bed', '.bim', '.fam'))
	file.symlink(paste0(input, i), paste0(output, i))

for (fail in fails) {
	write(sprintf('Working on %s ...', fail), stderr())
	if (isSampleFail(fail)) {
		params = list(
			bfile      = output,
			remove     = fail,
			out        = output,
			`make-bed` = T
		)
		cmd = sprintf("%s %s 1>&2", plink, cmdargs(params, equal = ' '))
		runcmd(cmd)
	}
}

for (f in Sys.glob(paste0(output, '.*~')))
	file.remove(f)