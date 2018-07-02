from pyppl import Box
from bioprocs.utils import runcmd, cmdargs, mem2

infile  = {{in.infile | repr}}
outfile = {{out.outfile | repr}}
umfile  = {{out.umfile | repr}}
tool    = {{args.tool | repr}}
picard  = {{args.picard | repr}}
chain   = {{args.lochain | repr}}
ref     = {{args.ref | repr}}
params  = {{args.params | repr}}
mem     = {{args.mem | repr}}
tmpdir  = {{args.tmpdir | repr}} 

# picard LiftoverVcf -Xmx4g -Xms1g  I=TCGA-05-4382-10.vcf O=1.vcf CHAIN=liftovers/hg38ToHg19.over.chain.gz R=ucsc_hg19.fa REJECT=r.vcf

if tool == 'picard':

	params.I      = infile
	params.O      = outfile
	params.CHAIN  = chain
	params.REJECT = umfile
	params.R      = ref

	javamem = mem2(mem, 'java')
	for jm in javamem.split()
		params[jm[1:]] = True

	params['Djava.io.tmpdir'] = tmpdir

	cmd = '{picard} LiftoverVcf {params}'
	runcmd(cmd.format(
		picard = picard,
		params = cmdargs(params, equal = '=', dash = '')
	))