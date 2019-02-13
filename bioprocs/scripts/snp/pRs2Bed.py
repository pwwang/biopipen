from os import path
from pyppl import Box
from bioprocs.utils import shell
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

snpfile   = {{ i.snpfile | quote}}
outfile   = {{ o.outfile | quote}}
inopts    = {{ args.inopts | repr}}
snpcol    = {{ args.snpcol or 0 | repr}}
dbsnp     = {{ args.dbsnp | quote}}
vcftools  = {{ args.vcftools | quote}}
sortby    = {{ args.sortby | quote}}
jobindir  = {{ job.indir | quote}}
joboutdir = {{ job.outdir | quote}}

# cruzdb
tool      = {{ args.tool | quote}}
genome    = {{ args.genome | quote}}
dbsnpver  = {{ args.dbsnpver | quote}}

if tool == 'cruzdb':
	'''	
	snp151(	
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

	# snps
	reader  = TsvReader(snpfile, cnames = False)
	snplist = list(set(r[snpcol] for r in reader))
	reader.close()

	from cruzdb import Genome
	g = Genome(genome)
	outfiletmp = outfile + '.tmp'
	writer = TsvWriter(outfiletmp)
	for i in range(0, len(snplist), 1000):
		chunk = snplist[i:i+1000]
		sql   = 'SELECT chrom, chromStart, chromEnd, name, score, strand, refUCSC, alleles, alleleFreqs FROM snp{dbsnpver} WHERE name in ({snps})'.format(
			dbsnpver = dbsnpver,
			snps     = ', '.join("'{}'".format(s) for s in chunk)
		)
		result = g.sql(sql)
		for r in result:
			allfreqs = dict(zip(r.alleles.split(','), r.alleleFreqs.split(',')))
			reffreq  = allfreqs.get(r.refUCSC, '0')
			if r.refUCSC in allfreqs:
				del allfreqs[r.refUCSC]
			if '' in allfreqs:
				del allfreqs['']
			writer.write([
				r.chrom, r.chromStart, r.chromEnd, r.name, r.score, r.strand, 
				r.refUCSC, ','.join(allfreqs.keys()), ','.join([reffreq] + list(allfreqs.values()))
			])
	writer.close()

else:
	# snps
	snplist = path.join(jobindir, path.basename(snpfile) + '.list')
	reader  = TsvReader(snpfile, cnames = False)
	writer  = TsvWriter(snplist)
	for r in reader:
		writer.write([r[snpcol]])
	reader.close()
	writer.close()

	shell.TOOLS.vcftools = vcftools
	vcftools = shell.Shell(equal = ' ', dash = '--').vcftools

	params = Box()
	params.snps = snplist
	params.recode = True
	params.out = path.join(joboutdir, 'tmp')
	if dbsnp.endswith('.gz'):
		params.gzvcf = dbsnp
	elif not path.isfile(dbsnp):
		raise ValueError('dbsnp file (args.dbsnp) is required by tool "local"')
	else:
		params.vcf = dbsnp
	vcftools(**params).run()

	reader = TsvReader(params.out + '.recode.vcf', cnames = False)
	outfiletmp = outfile + '.tmp'
	writer = TsvWriter(outfiletmp)
	for r in reader:
		writer.write([r[0], r[1], r[1], r[2], 0, '+', r[3], r[4]])
	reader.close()
	writer.close()

if sortby == 'coord':
	shell.sort(k = ['1,1', '2,2n'], _ = outfiletmp, _stdout = outfile)
else:
	shell.sort(k = '4', _ = outfiletmp, _stdout = outfile)

if tool == 'cruzdb':
	shell.rm_rf(outfiletmp)
else:
	shell.rm_rf(outfiletmp, params.out + '.recode.vcf')
