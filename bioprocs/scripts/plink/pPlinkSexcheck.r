library(methods)
{{rimport}}('__init__.r')

indir  = {{i.indir | R}}
outdir = {{o.outdir | R}}

plink    = {{args.plink | R}}

bedfile = Sys.glob(file.path(indir, '*.bed'))
input   = tools::file_path_sans_ext(bedfile)
output  = file.path(outdir, basename(input))

params = list(
	bfile       = input,
	`check-sex` = T,
	out         = output
)
cmd = sprintf("%s %s 1>&2", plink, cmdargs(params, equal = ' '))
runcmd(cmd)

sexcheck = read.table(paste0(output, '.sexcheck'), header = T, row.names = NULL, check.names = F)
sex.sample.fail = sexcheck[which(sexcheck$STATUS == 'PROBLEM'), c('FID', 'IID'), drop=F]
write.table(sex.sample.fail, paste0(output, '.sex.fail'), col.names = F, row.names = F, sep = "\t", quote = F)





