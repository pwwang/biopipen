"""Script for vcf.pVcfMutationAssessor"""

# pylint: disable=unused-import,invalid-name
from medoo import Medoo
from diot import Diot
from cyvcf2 import VCF, Writer

# pylint: disable=undefined-variable
infile = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
mutation_assessor_db = {{args.mutation_assessor_db | quote}}
# pylint: enable=undefined-variable

medoo = Medoo(dbtype='sqlite', database=f'file://{mutation_assessor_db}')

vcf = VCF(infile)
vcf.add_info_to_header({
    'ID': 'MutationAssessor',
    'Description': 'MutationAssessor annotations (gene|impact|score)',
    'Type': 'String',
    'Number': '1'
})

writer = Writer(outfile, vcf)

for record in vcf:
    # launch a query
    results = medoo.select(
        'scores',
        ['gene', 'impact', 'score'],
        where=dict(
            chrom=record.CHROM.replace('chr', ''),
            chrpos=record.POS,
            ref=record.REF,
            alt=record.ALT[0],
        )
    )
    if results.first():
        record.INFO['MutationAssessor'] = '|'.join(
            (results[0].gene, results[0].impact, str(results[0].score))
        )
    writer.write_record(record)

vcf.close()
writer.close()
