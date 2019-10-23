"""Commands of bcftools"""

from modkit import Modkit
from pyppl import Box, Proc
from . import params, delefactory, procfactory
Modkit().delegate(delefactory())

@procfactory
def _pQuery():
	"""
	@input
		infile: The input VCF file. Currently only single VCF file is supported.
	@output
		outfile: The output file
	@args:
		bcftools (str): Path to bcftools.
		params (Box): Other parameters for `bcftools query`.
			- See: https://samtools.github.io/bcftools/bcftools.html#query
	"""
	return Box(
		desc   = "Extracts fields from VCF or BCF files and outputs them in user-defined format",
		input  = 'infile:file',
		output = 'outfile:file:{{i.infile | stem | stem}}.query.txt',
		lang   = params.python.value,
		args   = Box(
			bcftools = params.bcftools.value,
			params = Box()
		)
	)

@procfactory
def _pView():
	"""
	@input
		infile: The input VCF file. Currently only single VCF file is supported.
		samfile: The sample name file to extract samples from infile
			- It also can be samples directly, separated by comma
			- If this is provided, `args.params.s` and `args.params.S` will be ignored
			- See https://samtools.github.io/bcftools/bcftools.html#view
	@output
		outfile: The output file
	@args:
		tabix    (str) : Path to tabix.
		bcftools (str) : Path to bcftools.
		gz       (bool): Whether output gzipped vcf file or not.
		nthread  (int) : Number of threads to use.
		params   (Box) : Other parameters for `bcftools view`.
			- See: https://samtools.github.io/bcftools/bcftools.html#view
	"""
	return Box(
		desc   = "View, subset and filter VCF or BCF files by position and filtering expression.",
		input  = 'infile:file, samfile:var',
		output = 'outfile:file:{{i.infile | stem | stem}}.vcf{{".gz" if args.gz else ""}}',
		lang   = params.python.value,
		args   = Box(
			gz       = False,
			nthread  = 1,
			tabix    = params.tabix.value,
			bcftools = params.bcftools.value,
			params   = Box()
		)
	)

@procfactory
def _pReheader():
	"""
	@input:
		infile: The input VCF file
		hfile: The new header file
		samfile: The new sample name file
			- It also can be samples directly, separated by comma
			- If this is provided, `args.params.samples` will be ignored
			- See https://samtools.github.io/bcftools/bcftools.html#reheader
	@output:
		outfile: The output VCF file with new header.
	@args:
		bcftools (str): Path to bcftools.
		nthread  (int): Number of threads to use.
		params   (Box): Other parameters for `bcftools view`.
			- See https://samtools.github.io/bcftools/bcftools.html#reheader
	"""
	return Box(
		desc   = "Modify header of VCF/BCF files, change sample names",
		input  = 'infile:file, hfile:file, samfile:var',
		output = 'outfile:file:{{i.infile | bn}}',
		lang   = params.python.value,
		args   = Box(
			bcftools = params.bcftools.value,
			nthread  = 1,
			params   = Box()
		)
	)

@procfactory
def _pFilter():
	"""
	@input:
		infile: The input VCF file
	@output:
		outfile: The output filtered VCF file
		statfile: Statistics of variants for each FILTER.
	@args:
		bcftools (str) : Path to bcftools
		nthread  (int) : Number of threads to use.
		gz       (bool): Whether output gzipped vcf file or not.
			- Overwrite `args.params.O`
		keep     (bool): Whether should we keep the excluded variants or not.
			- If not, args about filter names will not work.
		include  (str|list|dict): include only sites for which EXPRESSION is true.
			- See: https://samtools.github.io/bcftools/bcftools.html#expressions
			- Use `args.params.include` or `args.params.exclude` only if you just have one filter.
			- If provided, `args.params.include/exclude` will be ignored.
			- If `str`/`list` used, The filter names will be `Filter%d`
			- A dict is used when keys are filter names and values are expressions
		exclude  (str|list|dict): exclude sites for which EXPRESSION is true.
			- See also `args.include`
		params   (Box)          : Other parameters for `bcftools filter`
			- See: https://samtools.github.io/bcftools/bcftools.html#filter
	"""
	return Box(
		desc   = "Apply fixed-threshold filters to VCF files",
		input  = 'infile:file',
		output = 'outfile:file:{{i.infile | bn}}, statfile:file:{{i.infile | stem | stem}}.filterstats.txt',
		lang   = params.python.value,
		args   = Box(
			bcftools = params.bcftools.value,
			nthread  = 1,
			gz       = False,
			stat     = False,
			keep     = True,
			include  = None,
			exclude  = None,
			params   = Box(
				mode = '+' # accumulating Filter names instead of repolacing
			)
		)
	)

@procfactory
def _pAnnotate():
	"""
	@input:
		infile: The input VCF file
	@output:
		outfile: The annotated output VCF file
	@args:
		tabix   (str)  : Path to tabix, used to index `annfile`
		bcftools (str) : Path to bcftools
		annfile  (file): The annotation file.
			- See: https://samtools.github.io/bcftools/bcftools.html#annotate
		nthread (int)     : Number of threads to use.
		cols    (str|list): Overwrite `-c/--columns`.
		header  (str|list): headers to be added.
		params  (Box)     : Other parameters for `bcftools annotate`
	"""
	return Box(
		desc   = 'Add or remove annotations from VCF files',
		lang   = params.python.value,
		input  = 'infile:file',
		output = 'outfile:file:{{i.infile | bn}}',
		args   = Box(
			tabix    = params.tabix.value,
			bcftools = params.bcftools.value,
			nthread  = 1,
			annfile  = '',
			cols     = [],
			header   = [],
			params   = Box()
		)
	)
