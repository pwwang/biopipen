from pyppl import Box
from bioprocs.utils import mem2, logger, shell2 as shell

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
umfile  = {{o.umfile | quote}}
tool    = {{args.tool | quote}}
picard  = {{args.picard | quote}}
chain   = {{args.lochain | quote}}
ref     = {{args.ref | quote}}
params  = {{args.params | repr}}
mem     = {{args.mem | quote}}
tmpdir  = {{args.tmpdir | quote}} 

if not chain:
	logger.error('Chain file (args.lochain) not provided!')
	exit(1)

shell.load_config(picard = picard)

# picard LiftoverVcf -Xmx4g -Xms1g  I=TCGA-05-4382-10.vcf O=1.vcf CHAIN=liftovers/hg38ToHg19.over.chain.gz R=ucsc_hg19.fa REJECT=r.vcf

if tool == 'picard':

	params.I      = infile
	params.O      = outfile
	params.CHAIN  = chain
	params.REJECT = umfile
	params.R      = ref

	javamem = mem2(mem, 'java')
	for jm in javamem.split():
		params['-' + jm[1:]] = True

	params['-Djava.io.tmpdir'] = tmpdir
	shell.fg.picard.LiftoverVcf(**params)
