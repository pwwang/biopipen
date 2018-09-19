import sys
from cruzdb import Genome

# get the snps
snps = []
readopts = {{args.inopts}}
{% case True %}
{% when args.inmeta | isinstance: list %}
{{read.base.py | norepeats}}
reader = readBase({{i.infile | quote}}, **readopts)
reader.meta.add(*{{args.inmeta}})
{% when args.inmeta | isinstance: dict %}
{{read.base.py | norepeats}}
reader = readBase({{i.infile | quote}}, **readopts)
reader.meta.add(*{{args.inmeta}}.items())
{% when args.inmeta | bool %}
{{read, args.inmeta | lambda x, y: x[y]['py'] | norepeats}}
reader = locals()['read{{args.inmeta | @capitalize}}']({{i.infile | quote}}, **readopts)
{% else %}
ext = {{i.infile | ext | [1:] | quote}}
ext = ext[0].upper() + ext[1:]
reader = locals()['read' + ext]({{i.infile | quote}}, **readopts)
{% endif %}

for r in reader:
	r.START = int(r.START)
	r.END   = int(r.END)
	if r.START == r.END:
		r.END = r.START + 1
	snp = (r.CHR, r.START, r.END)
	if snp not in snps:
		snps.append(snp)

genome = {{args.genome | quote}}
g     = Genome (db=genome)
dbsnp = g.snp{{args.dbsnpver}}

{{write.bedx.py | norepeats}}
writer = writeBedx({{o.outfile | quote}})
writer.meta.add(*{{args.xcols}})
#writer.writeMeta()
writer.writeHead()
for snp in snps:
	chr, start, end = snp
	s = dbsnp.filter_by(chrom=chr, chromStart = start, chromEnd = end).first()
	'''
	snp147(
		bin=1566,
		chrom='chr8',
		chromStart=128700232L,
		chromEnd=128700233L,
		name='rs7005394',
		score=0,
		strand=u'+',
		refNCBI='T',
		refUCSC='T',
		observed='C/T',
		molType=u'genomic',
		class=u'single',
		valid=set(['by-2hit-2allele', 'by-cluster', 'by-hapmap', 'by-frequency', 'by-1000genomes']),
		avHet=0.498652,
		avHetSE=0.025923,
		func=set(['ncRNA']),
		locType=u'exact',
		weight=1L,
		exceptions=set([]),
		submitterCount=26,
		submitters='1000GENOMES,ABI,BCM-HGSC-SUB,BCM_SSAHASNP,BGI,BL,BUSHMAN,COMPLETE_GENOMICS,DDI,ENSEMBL,EVA-GONL,EVA_DECODE,EVA_GENOME_DK,EVA_UK10K_ALSPAC,EVA_UK10K_TWINSUK,GMI,HAMMER_LAB,HGSV,HUMANGENOME_JCVI,ILLUMINA-UK,JMKIDD_LAB,PJP,SSAHASNP,SSMP,TISHKOFF,WEILL_CORNELL_DGM,',
		alleleFreqCount=2,
		alleles='C,T,',
		alleleNs='2634.000000,2374.000000,',
		alleleFreqs='0.525958,0.474042,',
		bitfields=set(['maf-5-all-pops', 'maf-5-some-pop'])
	)
	'''
	if not s:
		sys.stderr.write('pyppl.log.warning: Cannot find snp at: %s\n' % str(snp))
	else:
		r        = readRecord()
		r.CHR    = s.chrom
		r.START  = str(s.chromStart)
		r.END    = str(s.chromEnd)
		r.NAME   = s.name
		r.SCORE  = str(s.score)
		r.STRAND = s.strand
		for xcol in {{args.xcols}}:
			setattr(r, xcol, str(getattr(s, xcol)))
		writer.write(r)
		