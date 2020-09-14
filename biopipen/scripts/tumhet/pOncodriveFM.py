"""Script for tumhet.pOncodriveFM"""

# pylint: disable=unused-import,invalid-name
import warnings
from pathlib import Path
from diot import Diot
from cyvcf2 import VCF
from bioprocs.utils import shell2 as shell

# pylint: disable=undefined-variable
infiles = {{i.infiles | repr}}
outgenes = {{o.outgenes | quote}}
outpathways = {{o.outpathways | quote}}
oncodrivefm = {{args.oncodrivefm | quote}}
gene2kegg = {{args.gene2kegg | quote}}
params = {{args.params | repr}}
# pylint: enable=undefined-variable

shell.load_config(oncodrivefm=oncodrivefm)

# compose input file, like:
# SAMPLE	GENE	SIFT	PPH2	MA
# 184	ENSG00000000971	0.81	0.973	0.46
# 280	ENSG00000001630	1	0.538	NA
# 30	ENSG00000001631	1	0.992	0.975

# gene => sample => annotator => score
input_data = dict()
annotators = set()

for infile in infiles:
    vcf = VCF(infile)

    if len(vcf.samples) != 1:
        warnings.warn(f"{infile} as more than 1 sample, "
                      "only the first one assumed.")

    sample = vcf.samples[0]
    for var in vcf:
        # We should get the same gene symbol from all annotations
        gene = None
        sifts = var.INFO.get('SIFTINFO')
        sift_score = None
        if sifts:
            annotators.add('SIFT')
            for sift in sifts.split(','):
                sift_items = sift.split('|')
                gene = gene or sift_items[3]
                if gene != sift_items[3]:
                    warnings.warn(f"{var} annotated with different "
                                  f"gene symbols: {gene}, {sift_items[3]}")
                score = sift_items[8]
                try:
                    score = float(score)
                except (TypeError, ValueError):
                    pass
                else:
                    sift_score = sift_score or score
                    sift_score = max(sift_score, score)
            if sift_score:
                input_data.setdefault(gene, {})
                input_data[gene].setdefault(sample, {})
                input_data[gene][sample]['SIFT'] = sift_score

        pp2 = var.INFO.get('PP2')
        if pp2:
            annotators.add('PPH2')
            pp2_items = pp2.split(',')
            gene = gene or pp2_items[0]
            pp2_score = pp2_items[4]
            try:
                pp2_score = float(pp2_score)
            except (TypeError, ValueError):
                pp2_score = None

            if pp2_score:
                input_data.setdefault(gene, {})
                input_data[gene].setdefault(sample, {})
                input_data[gene][sample]['PPH2'] = pp2_score

        ma = var.INFO.get('MutationAssessor')
        if ma:
            annotators.add('MA')
            ma_items = ma.split('|')
            gene = gene or ma_items[0]
            ma_score = ma_items[2]
            try:
                ma_score = float(ma_score)
            except (TypeError, ValueError):
                ma_score = None

            if ma_score:
                input_data.setdefault(gene, {})
                input_data[gene].setdefault(sample, {})
                input_data[gene][sample]['MA'] = ma_score

    vcf.close()

# gene => sample => annotator => score

# write input_data to tdm file
outgenes = Path(outgenes)
if not outgenes.stem.endswith('-genes'):
    raise ValueError('`o.outgenes` has to be named as `xxx-genes.tsv`')
tdmfile = outgenes.parent.joinpath(outgenes.stem[:-6] + '.tdm')
annotators = list(annotators)
with tdmfile.open('w') as ftdm:
    ftdm.write('\t'.join(['Sample', 'Gene'] + annotators) + '\n')
    for gene in input_data:
        # sample => annotator => score
        for sample, annoscore in input_data[gene].items():
            ftdm.write('\t'.join(
                [sample, gene] +
                [str(annoscore.get(anno, 'NA')) for anno in annotators]
            ) + '\n')

params.m = gene2kegg
params._ = tdmfile
params.o = outgenes.parent

shell.oncodrivefm(**params).fg()
