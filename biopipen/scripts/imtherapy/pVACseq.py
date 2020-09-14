from pathlib import Path
from diot import Diot
from bioprocs.utils import shell2 as shell
from bioprocs.utils.reference import vcfIndex
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter, TsvRecord

{% from os import path %}
{% from pyppl.utils import alwaysList %}
infile   = {{i.infile | quote}}
afile    = {{i.afile | ?path.isfile | =readlines | !alwaysList
					 | ?:len(_) == 1 and not _[0] | =:None | repr}}
outfile  = {{o.outfile | quote}}
outdir   = {{o.outdir | quote}}
bcftools = {{args.bcftools | quote}}
allele   = {{args.allele | ?path.isfile | =readlines | !alwaysList
						 | ?:len(_) == 1 and not _[0] | =:None | repr}}
netmhc     = {{args.netmhc | quote}}
iedb_mhc_i = {{args.iedb_mhc_i | quote}}
pvacseq    = {{args.pvacseq | quote}}
bdtool     = {{args.bdtool | quote}}
nthread    = {{args.nthread | quote}}
params     = {{args.params | repr}}

infile = vcfIndex(infile)

# get alleles
allele = afile or allele
if not allele:
	raise ValueError('No allele has been specified.')
allele = ','.join(allele)

shell.load_config(pvacseq = pvacseq, bcftools = bcftools)

bdtools = [	'MHCflurry','MHCnuggetsI','MHCnuggetsII','NNalign','NetMHC',
			'NetMHCIIpan','NetMHCcons','NetMHCpan','PickPocket','SMM','SMMPMBEC','SMMalign']
bdtools = {bdt.lower():bdt for bdt in bdtools}

# get sample name
sample = shell.bcftools.query(l = infile).splitlines()[0]

shell.rm_rf(Path(outdir).joinpath('MHC_Class_I', sample + '.tsv'), _debug = True)
params.t                      = nthread
params._                      = [infile, sample, allele, bdtools[bdtool], outdir]
params.k                      = params.get('k', True)
params.iedb_install_directory = Path(iedb_mhc_i).parent
shell.fg.pvacseq.run(**params)

# filter the epitopes with IC50(MT) >= 500 (SB) and IC50(WT) < 2000 (WB)
# Chromosome	Start	Stop	Reference	Variant	Transcript	Transcript Support Level	Ensembl Gene ID	Variant Type	Mutation	Protein Position	Gene Name	HGVSc	HGVSp	HLA Allele	Peptide Length	Sub-peptide Position	Mutation Position	MT Epitope Seq	WT Epitope Seq	Best MT Score Method	Best MT Score	Corresponding WT Score	Corresponding Fold Change	Tumor DNA Depth	Tumor DNA VAF	Tumor RNA Depth	Tumor RNA VAF	Normal Depth	Normal VAF	Gene Expression	Transcript Expression	Median MT Score	Median WT Score	Median Fold Change	NetMHC WT Score	NetMHC MT Score	cterm_7mer_gravy_score	max_7mer_gravy_score	difficult_n_terminal_residue	c_terminal_cysteine	c_terminal_proline	cysteine_count	n_terminal_asparagine	asparagine_proline_bond_count
reader = TsvReader(Path(outdir).joinpath('MHC_Class_I', sample + '.all_epitopes.tsv'))
writer = TsvWriter(outfile)
writer.cnames = ['HLA_allele', 'Peptide', 'Affinity', 'Gene', 'ENSG', 'ENST', 'Ref_peptide', 'Ref_affinity', 'Mutation', 'AAChange']
writer.writeHead()
for r in reader:
	out              = TsvRecord()
	out.HLA_allele   = r['HLA Allele']
	out.Peptide      = r['MT Epitope Seq']
	out.Affinity     = r['Best MT Score']
	out.Gene         = r['Gene Name']
	out.ENSG         = r['Ensembl Gene ID']
	out.ENST         = r['Transcript']
	out.Ref_peptide  = r['WT Epitope Seq']
	out.Ref_affinity = r['Corresponding WT Score']
	out.Mutation     = r.Chromosome + ':' + r.Start + '-' + r.Stop + '.' + r.Reference + '/' + r.Variant
	out.AAChange     = r.Mutation
	if float(out.Affinity) > 500 or float(out.Ref_affinity) < 2000:
		continue
	writer.write(out)
