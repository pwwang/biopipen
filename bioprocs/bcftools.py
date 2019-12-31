"""Commands of bcftools"""

from pyppl import Proc
from diot import Diot
from . import params, proc_factory

# @procfactory\ndef _(.+)\(\):\n\t("""[\s\S]+?""")\n\treturn Diot\(\n\t\t(desc\s+=.+)
# $1 = proc_factory(\n\t$3\n\tconfig = Diot(annotate = $2))

pQuery = proc_factory(
	desc   = "Extracts fields from VCF or BCF files and outputs them in user-defined format",
	lang   = params.python.value,
	config = Diot(annotate = """
	@input
		infile: The input VCF file. Currently only single VCF file is supported.
	@output
		outfile: The output file
	@args:
		bcftools (str): Path to bcftools.
		params (Diot): Other parameters for `bcftools query`.
			- See: https://samtools.github.io/bcftools/bcftools.html#query
	"""))
pQuery.input  = 'infile:file'
pQuery.output = 'outfile:file:{{i.infile | stem | stem}}.query.txt'
pQuery.args   = Diot(
	bcftools = params.bcftools.value,
	params   = Diot()
)

pView = proc_factory(
	desc   = "View, subset and filter VCF or BCF files by position and filtering expression.",
	config = Diot(annotate = """
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
		params   (Diot) : Other parameters for `bcftools view`.
			- See: https://samtools.github.io/bcftools/bcftools.html#view
	"""))
pView.input  = 'infile:file, samfile:var',
pView.output = 'outfile:file:{{i.infile | stem | stem}}.vcf{{".gz" if args.gz else ""}}',
pView.lang   = params.python.value,
pView.args   = Diot(
	gz       = False,
	nthread  = 1,
	tabix    = params.tabix.value,
	bcftools = params.bcftools.value,
	params   = Diot()
)

pReheader = proc_factory(
	desc   = "Modify header of VCF/BCF files, change sample names",
	config = Diot(annotate = """
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
		params   (Diot): Other parameters for `bcftools view`.
			- See https://samtools.github.io/bcftools/bcftools.html#reheader
	"""))
pReheader.input  = 'infile:file, hfile:file, samfile:var',
pReheader.output = 'outfile:file:{{i.infile | bn}}',
pReheader.lang   = params.python.value,
pReheader.args   = Diot(
	bcftools = params.bcftools.value,
	nthread  = 1,
	params   = Diot()
)

pFilter = proc_factory(
	desc   = "Apply fixed-threshold filters to VCF files",
	config = Diot(annotate = """
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
		params   (Diot)          : Other parameters for `bcftools filter`
			- See: https://samtools.github.io/bcftools/bcftools.html#filter
	"""))
pFilter.input  = 'infile:file',
pFilter.output = 'outfile:file:{{i.infile | bn}}, statfile:file:{{i.infile | stem | stem}}.filterstats.txt',
pFilter.lang   = params.python.value,
pFilter.args   = Diot(
	bcftools = params.bcftools.value,
	nthread  = 1,
	gz       = False,
	stat     = False,
	keep     = True,
	include  = None,
	exclude  = None,
	params   = Diot(
		mode = '+' # accumulating Filter names instead of repolacing
	)
)

pAnnotate = proc_factory(
	desc   = 'Add or remove annotations from VCF files',
	config = Diot(annotate = """
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
		params  (Diot)     : Other parameters for `bcftools annotate`
	"""))
pAnnotate.lang   = params.python.value,
pAnnotate.input  = 'infile:file',
pAnnotate.output = 'outfile:file:{{i.infile | bn}}',
pAnnotate.args   = Diot(
	tabix    = params.tabix.value,
	bcftools = params.bcftools.value,
	nthread  = 1,
	annfile  = '',
	cols     = [],
	header   = [],
	params   = Diot()
)

pConcat = proc_factory(
	desc   = 'Concatenate or combine VCF/BCF files with same samples in the same order.',
	config = Diot(annotate = """
	@input:
		infiles: The input vcf files
	@output:
		outfile: The output merged vcf file
	@args:
		nthread  (int) : The number of threads to use
		bcftools (path): The path to bcftools
		tabix    (path): The path to tabix, used to index vcf files.
		params   (Diot) : Other parameters for `bcftools concat`
		gz       (bool): Whether output gzipped vcf or not.
	"""))
pConcat.input  = 'infiles:files',
pConcat.output = 'outfile:file:{{i.infiles | [0] | stem2 | @append: "_etc.vcf"}}{{ \
	args.gz | ? | =:".gz" | !:""}}',
pConcat.lang   = params.python.value,
pConcat.args   = Diot(
	nthread  = 1,
	bcftools = params.bcftools.value,
	tabix    = params.tabix.value,
	params   = Diot(),
	gz       = False,
)
