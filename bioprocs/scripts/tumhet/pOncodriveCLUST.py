"""Script for tumhet.pOncodriveCLUST"""

# pylint: disable=unused-import,invalid-name
import re
import warnings
from pathlib import Path
from diot import Diot
from bioprocs.utils import shell2 as shell
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

# pylint: disable=undefined-variable
infile = {{i.infile | repr}}
outfile = {{o.outfile | quote}}
oncodriveclust = {{args.oncodriveclust | quote}}
oncodriveclust_data = {{args.oncodriveclust_data | quote}}
params = {{args.params | repr}}
# pylint: enable=undefined-variable

oncodrive_data = Path(oncodriveclust_data)
outfile = Path(outfile)

if not oncodrive_data.joinpath('gene_transcripts.tsv').exists():
    raise ValueError('`gene_transcripts.tsv` is required to be in '
                     '`oncodriveclust_data` directory')

shell.load_config(oncodriveclust=oncodriveclust)

# Make maf file to syn and nonsyn files for oncodriveclust
# nonsyn file should have mutation types: non-synonymous and stop
# syn file should only have synonymous mutations
MUTATION_CLASS_MAPS = dict(
    Silent='synonymous',
    Missense_Mutation='nonsynonymous',
    Nonsense_Mutation='stop',
    Nonstop_Mutation='nonsynonymous',
)

nonsyn_file = outfile.parent.joinpath('nonsyn.txt')
syn_file = outfile.parent.joinpath('syn.txt')

reader = TsvReader(infile)
nonsyn_writer = TsvWriter(nonsyn_file)
syn_writer = TsvWriter(syn_file)
nonsyn_writer.cnames = syn_writer.cnames = [
    'symbol', 'gene', 'transcript', 'Sample', 'ct', 'position'
]
nonsyn_writer.write_head()
syn_writer.write_head()

for r in reader:
    if r.Variant_Classification not in MUTATION_CLASS_MAPS:
        continue
    gene = r.Hugo_Symbol
    ensg = r.get('Gene',
                 r.get('HGNC_Ensembl_Gene_ID',
                       r.get('HGNC_Ensembl Gene ID'))).split('.')[0]
    enst = r.get('Transcript_ID', r.get('Annotation_Transcript')).split('.')[0]
    aachange = r.get('HGVSp_Short', r.get('Protein_Change'))
    sample = r.Tumor_Sample_Barcode
    matching = re.match(r'^p\.[^\d](\d+)[^\d]$', aachange)
    aapos = matching.group(1) if matching else ''

    if MUTATION_CLASS_MAPS[r.Variant_Classification] == 'synonymous':
        syn_writer.write([
            gene, ensg, enst, sample, 'synonymous', aapos
        ])
    else:
        nonsyn_writer.write([
            gene, ensg, enst, sample,
            MUTATION_CLASS_MAPS[r.Variant_Classification], aapos
        ])

reader.close()
nonsyn_writer.close()
syn_writer.close()

if oncodrive_data.joinpath('CGC_phenotype.tsv').exists():
    params.cgc = oncodrive_data.joinpath('CGC_phenotype.tsv')

if oncodrive_data.joinpath('pfam_domains.txt').exists():
    params.dom = oncodrive_data.joinpath('pfam_domains.txt')

params._ = [nonsyn_file, syn_file,
            oncodrive_data.joinpath('gene_transcripts.tsv')]
params.o = outfile

shell.oncodriveclust(**params).fg()
