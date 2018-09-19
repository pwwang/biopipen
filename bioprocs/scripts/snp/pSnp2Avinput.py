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
with open({{i.snpfile | quote}}) as f:
	snps = sorted(set([line.split(delimit)[col] for i, line in enumerate(f.read().splitlines()) if i >= skip and line.strip() and not line.startswith(comment)]))

genome = {{args.genome | quote}}

g     = Genome (db=genome)
dbsnp = g.snp{{args.dbsnpver}}

fout  = open ("{{o.outfile}}", "w")
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
		func    = list(s.func) # intron
		alleles = [a for a in s.alleles.split(',') if a and a!=ref] or ['N']

		#         chr pos pos ref all snp comm
		fout.write ("%s %s %s %s %s %s %s\n" % (
			s.chrom,
			s.chromEnd,
			s.chromEnd,
			ref,
			','.join(alleles),
			s.name,
			"/".join(func)
		))
fout.close()