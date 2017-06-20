from pyppl import proc

"""
Get SNP information from ucsc database
"""

"""
@name:
	pSNP2Bed
@input:
	`snpfile:file`: the snp file, each snp per line
@output:
	`outfile:file`: the result file, columns are:
	- chrom, start(0-based), end, name, score, strand, ref, allele
@args:
	`genome`: default: hg19
	`snpver`: default: snp147
@requires:
	[`python-cruzdb`](https://github.com/brentp/cruzdb)
"""
pSNP2Bed = proc()
pSNP2Bed.input  = "snpfile:file"
pSNP2Bed.output = "outfile:file:{{snpfile | fn}}.bed"
pSNP2Bed.args   = {"genome": "hg19", "snpver": "snp147"}
pSNP2Bed.lang   = "python"
pSNP2Bed.script = """
from cruzdb import Genome
import sys

g     = Genome (db="{{proc.args.genome}}")
dbsnp = g.{{proc.args.snpver}}
snps  = list(set([line.split()[0] for line in open("{{snpfile}}") if line.strip()]))
fout  = open ("{{outfile}}", "w")
for snp in snps:
	s = dbsnp.filter_by(name=snp).first()
	if not s:
		sys.stderr.write(s + "\\n")
	else:
		#           chr  st   end  name score strd ref allele
		ns = s.alleleNs.split(',')
		ns = filter(None, ns)
		ns = map(float, ns)
		i  = ns.index(max(ns))
		al = s.alleles.split(',')
		al = filter(None, al)
		a  = al[i]
		r  = s.refNCBI
		if a == r and len(ns) > 1:
			i = len(ns) - 1 - i
			a = al[i]
		func = map (str, list(s.func))
		fout.write ("%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" % (
			s.chrom,
			s.chromStart,
			s.chromEnd,
			s.name,
			s.score,
			s.strand,
			r,
			a,
			'/'.join(func)
		))
fout.close()
"""

"""
@name:
	pSNP2Avinput
@input:
	`snpfile:file`: the snp file, each snp per line
@output:
	`outfile:file`: the result avinput file
@args:
	`genome`: default: hg19
	`snpver`: default: snp147
@requires:
	[`python-cruzdb`](https://github.com/brentp/cruzdb)
"""
pSNP2Avinput = proc ()
pSNP2Avinput.input  = "snpfile:file"
pSNP2Avinput.output = "outfile:file:{{snpfile | fn}}.avinput"
pSNP2Avinput.lang   = "python"
pSNP2Avinput.args   = {"genome": "hg19", "snpver": "snp147"}
pSNP2Avinput.script = """
from cruzdb import Genome
import sys

g     = Genome (db="{{proc.args.genome}}")
dbsnp = g.{{proc.args.snpver}}
snps  = list(set([line.split()[0] for line in open("{{snpfile}}") if line.strip()]))
fout  = open ("{{outfile}}", "w")
for snp in snps:
	s = dbsnp.filter_by(name=snp).first()
	if not s:
		sys.stderr.write(s + "\\n")
	else:
		ns = s.alleleNs.split(',')
		ns = filter(None, ns)
		ns = map(float, ns)
		i  = ns.index(max(ns))
		al = s.alleles.split(',')
		al = filter(None, al)
		a  = al[i]
		r  = s.refNCBI
		if a == r and len(ns) > 1:
			i = len(ns) - 1 - i
			a = al[i]
		func = map (str, list(s.func))

		#         chr pos pos ref all snp comm
		fout.write ("%s %s %s %s %s %s %s\\n" % (
			s.chrom,
			s.chromEnd,
			s.chromEnd,
			r,
			a,
			s.name,
			"/".join(func)
		))
fout.close()
"""

