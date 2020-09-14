from pathlib import Path
from diot import Diot
from bioprocs.utils import shell2 as shell
from bioprocs.utils.reference import bamIndex
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile   = {{i.infile | quote}}
outdir   = Path({{o.outdir | quote}})
outfile  = {{o.outfile | quote}}
sample   = {{i.infile | stem2 | quote}}
hla_la   = {{args.hla_la | quote}}
picard   = {{args.picard |quote}}
bwa      = {{args.bwa | quote}}
java     = {{args.java | quote}}
samtools = {{args.samtools | quote}}
nthread  = {{args.nthread | repr}}
params   = {{args.params | repr}}

bamIndex(infile, samtools = samtools)
shell.load_config(hla_la = hla_la)

params.picard_sam2fastq_bin = shell.which(picard)
params.bwa_bin              = shell.which(bwa)
params.samtools_bin         = shell.which(samtools)
params.java_bin             = shell.which(java)
params.maxThreads           = nthread
params.workingDir           = outdir.parent
params.sampleID             = sample
params.BAM                  = infile
params.graph                = 'PRG_MHC_GRCh38_withIMGT'

shell.hla_la(**params).fg

hlafile = outdir.joinpath('hla', 'R1_bestguess_G.txt')
reader = TsvReader(hlafile)
writer = TsvWriter(outfile)
writer.cnames = ['Type', 'Allele', 'Reported'] + [cname for cname in reader.cnames if cname != 'Allele']
writer.writeHead()

hits = []
for r in reader:
	if r.Locus not in hits:
		hits.append(r.Locus)
		r.Type = r.Locus + ('_' if len(r.Locus) > 1 else '') + '1'
	else:
		r.Type = r.Locus + ('_' if len(r.Locus) > 1 else '') + '2'
	r.Reported = r.Allele
	r.Allele = 'HLA-' + r.Allele
	writer.write(r)
