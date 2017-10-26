import sys
from cruzdb import Genome

{% if args.header %}
skip = {{args.skip}} + 1
{% else %}
skip = {{args.skip}}
{% endif %}

# get the snps
delimit = {{args.delimit | quote}}
comment = {{args.comment | quote}}
col     = {{args.col}}
with open({{in.snpfile | quote}}) as f:
	snps = sorted(set([line.split(delimit)[col] for i, line in enumerate(f.read().splitlines()) if i >= skip and line.strip() and not line.startswith(comment)]))

genome = {{args.genome | quote}}

g     = Genome (db=genome)
dbsnp = g.snp{{args.dbsnpver}}

fout  = open ("{{out.outfile}}", "w")
for snp in snps:
	s = dbsnp.filter_by(name=snp).first()
	if not s:
		sys.stderr.write('pyppl.log.warning: Cannot find coordinates for SNP: %s\n' % snp)
	else:
		# chr start end  name score strand otherinfo
		chrom   = s.chrom
		start   = s.chromStart
		end     = s.chromEnd
		name    = snp
		strand  = s.strand
		ref     = s.refUCSC
		avHet   = s.avHet # average Het
		avHetSE = s.avHetSE # average Het standard error
		func    = list(s.func) # intron
		alleles = s.alleles.split(',')
		allfrqs = s.alleleFreqs.split(',')
		allouts = '|'.join([a + '$' + str(allfrqs[i]) for i,a in enumerate(alleles) if a and a!=ref])

		fout.write ("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
			chrom,
			start,
			end,
			name,
			0,
			strand,
			'; '.join([
				'REF:%s' % ref,
				'AVHET:%s' % avHet,
				'AVHETSE:%s' % avHetSE,
				'FUNC:%s' % ','.join(func),
				'ALLELES:%s' % allouts
			])
		))
fout.close()