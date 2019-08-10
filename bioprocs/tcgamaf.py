"""Operations on TCGA MAF file"""
from pyppl import Box
from . import params, delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pMaf2Vcf():
	"""
	@input:
		infile: The input MAF file
	@output:
		outfile: Output multi-sample VCF containing all TN-pairs
		outdir: A directory with output VCF files
	@args:
		maf2vcf (str): Path to `maf2vcf.pl` from `vcf2maf` tool.
		merge: Whether merge the VCF for individual samples or not.
		params: Other parameters for `maf2vcf.pl`
	"""
	return Box(
		desc = 'Convert MAF file back to VCF files',
		lang = params.python.value,
		input = 'infile:file',
		output = [
			'outfile:file:{{i.infile | stem}}.vcfs/{{i.infile | stem}}.vcf',
			'outdir:dir:{{i.infile | stem}}.vcfs'],
		args = Box(
			maf2vcf = params.maf2vcf.value,
			ref     = params.ref.value,
			params  = Box({'per-tn-vcfs': True})),
	)


@procfactory
def _pMafAddChr():
	"""
	@input:
		infile: The input MAF file
	@output:
		outfile: The output MAF file
	"""
	return Box(
		desc = 'Add chr to chromosome if not present',
		lang = params.python.value,
		input = 'infile:file',
		output = 'outfile:file:{{i.infile | stem}}.maf'
	)