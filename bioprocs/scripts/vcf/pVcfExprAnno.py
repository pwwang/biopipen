from pathlib import Path
from diot import Diot
from bioprocs.utils import shell2 as shell
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile                   = {{i.infile | quote}}
exprfile                 = {{i.exprfile | quote}}
outfile                  = Path({{o.outfile | quote}})
vcf_expression_annotator = {{args.vcf_expression_annotator | quote}}
exprtype                 = {{args.exprtype | quote}}
params                   = {{args.params | repr}}
istx                     = {{args.istx | repr}}
bcftools                 = {{args.bcftools | quote}}
sample                   = {{args.sample | repr}}

shell.load_config(vcf_expression_annotator = vcf_expression_annotator, bcftools = bcftools)
if isinstance(sample, str):
	params.s = sample
else:
	samples = shell.bcftools.query(l = infile).splitlines()
	samples = [line.strip() for line in samples if line.strip()]
	if len(samples)>1:
		params.s = samples[sample]

if exprtype == 'kallisto':
	reader         = TsvReader(exprfile)
	efile_withgene = outfile.parent.joinpath(outfile.stem + '.gene.expr.txt')
	writer         = TsvWriter(efile_withgene)
	writer.cnames  = reader.cnames[:-1] + ['abundance', 'gene_name']
	writer.writeHead()
	for r in reader:
		gene = r[0] if '::' not in r[0] else r[0].split('::', 1)[0]
		r.abundance = r.tpm
		r.gene_name = gene
		writer.write(r)
	writer.close()
	exprfile = efile_withgene

params.o = outfile
shell.vcf_expression_annotator(
	params, infile, exprfile, exprtype, 'transcript' if istx else 'gene'
).fg
