from bioprocs.utils import runcmd, cmdargs
from bioprocs.utils.tsvio import TsvReader, TsvWriter, TsvRecord

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
liftover = {{args.liftover | quote}}
lochain  = {{args.lochain | quote}}
genome   = {{args.genome | quote}} # target genome

bedfile      = "{{job.outdir}}/{{i.infile | fn}}.bed"
mappedfile   = "{{job.outdir}}/{{i.infile | fn}}.mapped.bed"
unmappedfile = "{{job.outdir}}/{{i.infile | fn}}.unmapped.bed"

# get a bed file from maf file
reader = TsvReader(infile, ftype = 'head')
writer = TsvWriter(bedfile, ftype = 'bed')

for i, r in enumerate(reader):
	bedr        = TsvRecord()
	bedr.CHR    = r.Chromosome
	bedr.START  = r.Start_Position
	bedr.NAME   = 'MAF{}'.format(i)
	bedr.SCORE  = int(r.End_Position) - int(r.Start_Position)
	bedr.STRAND = r.Strand
	# liftOver doesn't work with bed index start with 0
	bedr.END = int(r.End_Position) + 1 if bedr.SCORE == 0 else r.End_Position
	writer.write(bedr)

writer.close()
reader.rewind()

# run liftOver
cmd = '{liftover} {bedfile} {chainfile} {mappedfile} {unmappedfile}'
runcmd(cmd.format(
	liftover     = liftover,
	bedfile      = bedfile,
	chainfile    = lochain,
	mappedfile   = mappedfile,
	unmappedfile = unmappedfile
))

# read all mapped file
reader2 = TsvReader(mappedfile, ftype = 'bed')
mapped = {}
for r in reader2:
	if r.SCORE == 0:
		mapped[r.NAME] = r.CHR, r.START, r.START
	else:
		mapped[r.NAME] = r.CHR, r.START, r.END

# mapped them back
writer2 = TsvWriter(outfile)
writer2.meta = reader.meta
writer2.writeHead()
for i, r in enumerate(reader):
	name = 'MAF{}'.format(i)
	if not name in mapped: continue
	Chr, Start, End = mapped[name]
	r.Chromosome, r.Start_Position, r.End_Position = mapped[name]
	r.NCBI_Build = genome
	writer2.write(r)
		
