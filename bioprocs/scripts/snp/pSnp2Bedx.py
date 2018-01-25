import sys
from cruzdb import Genome

# get the snps
snps = []
readopts = {{args.inopts}}
{% if args.inmeta | lambda x: isinstance(x, list) %}
{{read.base.py}}
reader = readBase({{in.snpfile | quote}}, **readopts)
reader.meta.add(*{{args.inmeta}})
{% elif args.inmeta | lambda x: isinstance(x, dict) %}
{{read.base.py}}
reader = readBase({{in.snpfile | quote}}, **readopts)
reader.meta.add(*{{args.inmeta}}.items())
{% elif args.inmeta %}
{{read, args.inmeta | lambda x, y: x[y]['py']}}
reader = locals()['read{{args.inmeta | lambda x: x[0].upper() + x[1:]}}']({{in.snpfile | quote}}, **readopts)
{% else %}
ext = {{in.snpfile | ext | [1:] | quote}}
ext = ext[0].upper() + ext[1:]
reader = locals()['read' + ext]({{in.snpfile | quote}}, **readopts)
{% endif %}

for r in reader:
	if r.SNP not in snps:
		snps.append(r.SNP)

genome = {{args.genome | quote}}
g     = Genome (db=genome)
dbsnp = g.snp{{args.dbsnpver}}

{{write.bedx.py}}
writer = writeBedx({{out.outfile | quote}})
writer.meta.add(*[x.upper() for x in {{args.xcols}}])
writer.writeMeta()
writer.writeHead()
for snp in snps:
	s = dbsnp.filter_by(name=snp).first()
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
		sys.stderr.write('pyppl.log.warning: Cannot find coordinates for SNP: %s\n' % snp)
	else:
		r        = readRecord()
		r.CHR    = s.chrom
		r.START  = str(s.chromStart)
		r.END    = str(s.chromEnd)
		r.NAME   = snp
		r.SCORE  = str(s.score)
		r.STRAND = s.strand
		for xcol in {{args.xcols}}:
			setattr(r, xcol.upper(), str(getattr(s, xcol)))
		writer.write(r)