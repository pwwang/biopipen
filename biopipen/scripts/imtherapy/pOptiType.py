from pathlib import Path
from diot import Diot
from bioprocs.utils import shell2 as shell
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

fqfile1  = Path({{i.fqfile1 | quote}})
fqfile2  = {{i.fqfile2 | quote}} # not sure we have it or not
outdir   = Path({{o.outdir | quote}})
optitype = {{args.optitype | quote}}
picard   = {{args.picard | quote}}
bwa      = {{args.bwa | quote}}
hlaref   = {{args.hlaref | quote}}
params   = {{args.params | repr}}
nthread  = {{args.nthread | repr}}

bwaps = Diot(params.pop('bwa', {}))

shell.load_config(picard = picard, bwa = bwa, optitype = optitype)

# Extract HLA sequences
# Let's check if we have a bam file or a pair of pair-end fastq file
if fqfile1.suffix == '.bam':
	# convert it to fastq
	tmpfq1 = outdir.joinpath(outdir.stem + '_1.fq')
	tmpfq2 = outdir.joinpath(outdir.stem + '_2.fq')
	opts = Diot(F = tmpfq1, F2 = tmpfq2, I = fqfile1)
	shell.picard.SamToFastq(**opts).fg
	fqfile1 = tmpfq1
	fqfile2 = tmpfq2

# Use bwa to map it to hla reference
tmpbam = outdir.joinpath(outdir.stem + '.sam')

bwaps.o = tmpbam
bwaps._ = [hlaref, fqfile1, fqfile2]
bwaps.t = nthread
shell.bwa.mem(**bwaps).fg

fqfish1 = outdir.joinpath(outdir.stem + '_fished_1.fq')
fqfish2 = outdir.joinpath(outdir.stem + '_fished_2.fq')
opts = Diot(I = tmpbam, F = fqfish1, F2 = fqfish2)
shell.picard.SamToFastq(**opts).fg

params.i = [fqfish1, fqfish2]
params.verbose = True
params.outdir = outdir
shell.optitype(**params).fg

# replace this:
# 	A1	A2	B1	B2	C1	C2	Reads	Objective
# 0	A*24:02	A*02:01	B*44:02	B*35:01	C*12:03	C*07:27	444.0	431.99200000000184
# to:
# Type	Allele	Reads	Objective
# A1	HLA-A*24:02	444.0	431.99
# A2	HLA-A*02:01	444.0	431.99
# ...

sample        = outdir.stem
resfile       = list(outdir.glob('*/*_result.tsv'))[0]
reader        = TsvReader(resfile)
writer        = TsvWriter(outdir.joinpath(sample + '.optitype.txt'))
writer.cnames = ['Type', 'Allele', 'Reads', 'Objective']
writer.writeHead()
for r in reader:
	for subtype in ('A1', 'A2', 'B1', 'B2', 'C1', 'C2'):
		writer.write([subtype, 'HLA-' + r[subtype], r.Reads, r.Objective])
