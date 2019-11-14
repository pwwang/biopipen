"""A set of processes to generate/process vcf files"""
from os import path
from glob import glob
from pyppl import Proc, Box
from . import params
from .utils import fs2name
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pVcfFilter():
	"""
	@name:
		pVcfFilter
	@description:
		Filter records in vcf file.
	@input:
		`infile:file`: The input file
	@output:
		`outfile:file`: The output file
	@args:
		`filters`: A dict of filters like `{<filtername>: <filter>}`. `<filter>` should be a string of lambda function:
			```
			"lambda record, samples: <expression>"
			* ``record.CHROM`` : 'chr20'
			* ``record.POS``   : 1234567
			* ``record.ID``    : 'microsat1'
			* ``record.REF``   : ''GTC''
			* ``record.ALT``   : [G, GTCT]
			* ``record.QUAL``  : 50
			* ``record.FILTER``: ['PASS'] # NO!, PASS should be []
			* ``record.INFO``  : {'AA': 'G', 'NS': 3, 'DP': 9}
			* samples = record.samples
			* len(samples): 3
			* samples[0].sample: 'NA00001'
			* samples[0]: Call(sample=NA00001, CallData(GT=0/1, GQ=35, DP=4))
			* samples[0].data: calldata(GT='0/1', GQ=35, DP=4)
			* samples[0].data.GT: '0/1'
			```
			- see here for record and samples: https://github.com/jamescasbon/PyVCF
			- Remember if filters() returns True, record filtered.
			- For builtin filters, you may specify them as `{<filter>: <param>}`
			- You can also use `!` to specify a negative builtin filter: `{!<filter>: <param>}`
			- Bulitin filters:
				- SNPONLY: keeps only SNPs (`{"!SNPONLY": None}` means filtering SNPs out)
				- BIALTONLY: keeps only mutations with bi-allele
				- QUAL: keeps only site quality >= param (`{'QUAL': 30}`)
		`gz`     : Whether to gzip the output file. Default: False
		`keep`   : Whether to keep the filtered records. Default: True. (only for gatk, snpsift at filter step)
	@requires:
		[`pyvcf`](https://github.com/jamescasbon/PyVCF)
	"""
	pVcfFilter              = Proc(desc = 'Filter records in vcf file.')
	pVcfFilter.input        = "infile:file"
	pVcfFilter.output       = "outfile:file:{{i.infile | fn2}}.vcf{% if args.gz %}.gz{% endif %}"
	pVcfFilter.args.filters = Box()
	pVcfFilter.args.gz      = False
	pVcfFilter.args.keep    = True # only for gatk, snpsift at filter step
	pVcfFilter.lang         = params.python.value
	pVcfFilter.script       = "file:scripts/vcf/pVcfFilter.py"
	return pVcfFilter

@procfactory
def _pVcfUnique():
	"""
	@name:
		pVcfUnique
	@description:
		Remove duplicate mutations from a VCF file.
		Because in most case where we want to remove the duplicate mutations, it might be
		because other program not accepting them. In this case, we don't put a filter on
		the records, but just remove them instead.
	@input:
		`infile:file`: The input vcf file.
	@output:
		`outfile:file`: The output vcf file. Default: `{{i.infile | fn2}}.vcf{% if args.gz %}.gz{% endif %}`
	@args:
		`upart`: The unique parts. Could be part of: `['CHROM', 'POS', 'ID', 'REF', 'ALT']`
		`keep` : Which record to keep.
			- `bisnp` : Snp with bi-allele
			- `snp`   : Snps
			- `bialt` : Bi-allele mutations
			- `first` : The first record (default)
			- `last`  : The last record
			- `random`: A random record
			- Multiple ways can be used: `first, snp` is to select first snp (in case multiple snps duplicated)
		`gz`: Bgzip the output vcf file or not. Default: `False`
	@requires:
		`pyvcf`
	"""
	pVcfUnique            = Proc(desc = 'Remove duplicate mutations from a VCF file.')
	pVcfUnique.input      = 'infile:file'
	pVcfUnique.output     = 'outfile:file:{{i.infile | fn2}}.vcf{% if args.gz %}.gz{% endif %}'
	pVcfUnique.args.upart = ['CHROM', 'POS', 'ID', 'REF', 'ALT']
	pVcfUnique.args.keep  = 'first' # last, random, snp, bisnp, bialt
	pVcfUnique.args.gz    = False
	pVcfUnique.lang       = params.python.value
	pVcfUnique.script     = "file:scripts/vcf/pVcfUnique.py"
	return pVcfUnique

@procfactory
def _pVcfRemoveFilter():
	"""
	@name:
		pVcfRemoveFilter
	@description:
		Remove one or more filters in vcf files
	@input:
		`infile:file`: The input vcf file
	@output:
		`outfile:file`: The output file
	@args:
		`rmfilter`: The filters to remove. If None, ALL filters will be removed!
			- A `list` of filter names.
	"""
	pVcfRemoveFilter               = Proc(desc = 'Remove one or more filters in vcf files')
	pVcfRemoveFilter.input         = 'infile:file'
	pVcfRemoveFilter.output        = 'outfile:file:{{i.infile | bn}}'
	pVcfRemoveFilter.args.rmfilter = None # remove all filters
	pVcfRemoveFilter.lang          = params.python.value
	pVcfRemoveFilter.script        = "file:scripts/vcf/pVcfRemoveFilter.py"
	return pVcfRemoveFilter

@procfactory
def _pVcf():
	"""
	@name:
		pVcf
	@description:
		Use pyvcf to manipulate vcf file
	@input:
		`infile:file`: The input vcf file
	@output:
		`outfile:file`: The output vcf file
	@args:
		`helper`: The helper code injected to script
			- Since lambda function can't do assignment and manipulation so you can write some help function here
		`reader`: A string of lambda function to manipulate the reader (arg: vcf.Reader instance)
		`record`: A string of lambda function to manipulate the record (arg1: vcf.Record instance, arg2: vcf.Writer, arg3: `samples`)
		`gz`: Gzip the ouput file
	"""
	pVcf             = Proc(desc = 'Motify vcf file using pyvcf')
	pVcf.input       = "infile:file"
	pVcf.output      = "outfile:file:{{i.infile | fn}}.vcf{% if args.gz %}.gz{% endif %}"
	pVcf.args.helper = ''
	pVcf.args.reader = None
	pVcf.args.record = None
	pVcf.args.gz     = False
	pVcf.lang        = params.python.value
	pVcf.script      = "file:scripts/vcf/pVcf.py"
	return pVcf

@procfactory
def _pVcfAnno():
	"""
	@name:
		pVcfAnno
	@description:
		Annotate the variants in vcf file.
		You have to prepare the databases for each tool.
	@input:
		`infile:file`: The input vcf file
	@output:
		`outfile:file`: The output file (output file of annovar will also be converted to vcf)
		`outdir`: The output directory, used to fetch some stat/summary files
	@args:
		`tool`:            The tool used to do annotation. Default: snpeff
		`snpeff`:          The path of snpeff. Default: snpEff
		`vep`:             The path to vep. Default: vep
		`gz`:              Whether to gzip the result file. Default: False
		`annovar`:         The path of annovar. Default: annotate_variation.pl
		`annovar_convert`: The path of convert2annovar.pl, used to convert vcf to annovar input file. Default: convert2annovar.pl
		`genome`:          The genome for annotation. Default: hg19
		`tmpdir`:          The tmpdir, mainly used by snpeff. Default: <system tmpdir>
		`dbs`:             The path of database for each tool. Required by 'annovar' and 'vep'
		`params`:          Other params for tool. Default: ''
		`snpeffStats`:     Whether to generate stats file when use snpeff. Default: False
		`mem`:             The memory used by snpeff. Default: '4G'
	@requires:
		[`annovar`](http://doc-openbio.readthedocs.io/projects/annovar/en/latest/)
		[`snpeff`](http://snpeff.sourceforge.net/SnpEff_manual.html#intro)
		[`vep`](http://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html)
	"""
	pVcfAnno                      = Proc(desc = 'Annotate the variants in vcf file.')
	pVcfAnno.input                = "infile:file"
	pVcfAnno.output               = [
		"outfile:file:{{i.infile | fn}}.{{args.tool}}/{{i.infile | fn}}.{{args.tool}}.vcf{% if args.gz %}.gz{% endif %}",
		"outdir:dir:{{i.infile | fn}}.{{args.tool}}"
	]
	pVcfAnno.echo                 = Box(jobs = 0, type = 'stderr')
	pVcfAnno.args.tool            = 'snpeff'
	pVcfAnno.args.snpeff          = params.snpeff.value
	pVcfAnno.args.vep             = params.vep.value
	pVcfAnno.args.gz              = False
	pVcfAnno.args.vcfanno         = params.vcfanno.value
	pVcfAnno.args.annovar         = params.annovar.value
	pVcfAnno.args.annovar_convert = params.annovar_convert.value
	pVcfAnno.args.genome          = params.genome.value
	pVcfAnno.args.tmpdir          = params.tmpdir.value
	pVcfAnno.args.dbs             = Box(
		snpeff  = params.snpeffDb.value,
		annovar = params.annovarDb.value,
		vep     = params.vepDb.value,
		vcfanno = []
	)
	pVcfAnno.args.snpeffStats = False
	pVcfAnno.args.nthread     = 1
	pVcfAnno.args.params      = Box()
	pVcfAnno.args.mem         = params.mem8G.value
	pVcfAnno.lang             = params.python.value
	pVcfAnno.script           = "file:scripts/vcf/pVcfAnno.py"
	return pVcfAnno

@procfactory
def _pVcfSplit():
	"""
	@name:
		pVcfSplit
	@description:
		Split multi-sample Vcf to single-sample Vcf files.
	@input:
		`infile:file`: The input vcf file
		`samples`:     The samples, if not provided, will extract all samples
	@output:
		`outdir:dir`:  The output directory containing the extracted vcfs
	@args:
		`tool`:     The tool used to do extraction. Default: bcftools (gatk, awk)
		`bcftools`: The path of bcftools, used to extract the sample names from input vcf file.
		`gatk`:     The path of gatk.
	"""
	pVcfSplit                     = Proc(desc = "Split multi-sample Vcf to single-sample Vcf files.")
	pVcfSplit.input               = "infile:file, samples"
	pVcfSplit.output              = "outdir:dir:{{i.infile | fn}}-individuals"
	pVcfSplit.args.tool           = 'bcftools'
	pVcfSplit.args.bcftools       = params.bcftools.value # used to extract samples
	pVcfSplit.args.gatk           = params.gatk.value
	pVcfSplit.args.tabix          = params.tabix.value
	pVcfSplit.args.ref            = params.ref.value # only for gatk
	pVcfSplit.args.params         = Box()
	pVcfSplit.args.nthread        = 1
	pVcfSplit.lang                = params.python.value
	pVcfSplit.script              = "file:scripts/vcf/pVcfSplit.py"
	return pVcfSplit

@procfactory
def _pVcfMerge():
	"""
	@name:
		pVcfMerge
	@description:
		Merge single-sample Vcf files to multi-sample Vcf file.
	@input:
		`infiles:files`: The input vcf files
		`outfile:dir`:  The output multi-sample vcf.
	@args:
		`tool`:     The tool used to do extraction. Default: bcftools
		`vcftools`: The path of vcftools' vcf-subset
		`bcftools`: The path of bcftools, used to extract the sample names from input vcf file.
		`gatk`:     The path of gatk.
	"""
	pVcfMerge               = Proc(desc = "Merge single-sample Vcf files to multi-sample Vcf file.")
	pVcfMerge.input         = "infiles:files"
	pVcfMerge.output        = "outfile:file:{{i.infiles | fs2name}}.vcf{{'.gz' if args.gz else ''}}"
	pVcfMerge.args.tool     = 'bcftools'
	pVcfMerge.args.vcftools = params.vcftools_merge.value
	pVcfMerge.args.bcftools = params.bcftools.value
	pVcfMerge.args.gatk     = params.gatk.value
	pVcfMerge.args.params   = Box()
	pVcfMerge.args.tabix    = params.tabix.value
	pVcfMerge.args.ref      = params.ref.value # only for gatk
	pVcfMerge.args.gz       = False
	pVcfMerge.args.nthread  = 1
	pVcfMerge.envs.fs2name  = fs2name
	pVcfMerge.lang          = params.python.value
	pVcfMerge.script        = "file:scripts/vcf/pVcfMerge.py"
	return pVcfMerge

@procfactory
def _pVcf2Maf():
	"""
	@input:
		`infile:file` : The input vcf file
			- see `args.tumor`
	@output:
		`outfile:file`: The output maf file
	@args:
		`tool`     : Which tool to use, vcf2maf or oncotator.
			- Only hg19 available for oncotator
		`vcf2maf`  : The path of vcf2maf.pl
		`vep`      : The path of vep
			- For `vcf2maf`
		`vepDb`    : The path of database for vep
		`filtervcf`: The filter vcf. Something like: ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
			- For `vcf2maf`
		`ref`      : The reference genome
		`nthread`  : Number of threads used to extract samples.
		`bcftools` : Path to bcftools used to extract sample names.
		`tumor`    : The index of the tumor sample in the vcf file. Default: `auto`
			- `auto`: The sample name matches the file name as the tumor, only if there are 2 samples.
			- `0`: The first sample as the tumor.
			- `1`: The second sample as the tumor.
			- Otherwise: All samples treated as tumors.
			- If there is only one sample, this option has no effect
		tabix: Path to tabix, used to index Vcf file.
		oncotator: Path to oncotator.
		oncotator_db (dir): Path to oncotator database.
		params (Box): Extra parameters for the tool.
		withchr (bool): Should we add chr to Chromosome column or not.
		genome (str): The genome used to replace __UNKNOWN__ for the NCBI_Build column
	"""
	return Box(
		desc   = 'Convert Vcf file to Maf file',
		input  = 'infile:file',
		output = 'outfile:file:{{i.infile | fn2}}.maf',
		lang   = params.python.value,
		args   = Box(
			tool         = 'oncotator',
			withchr      = True,
			vcf2maf      = params.vcf2maf.value,
			vep          = params.vep.value,
			tabix        = params.tabix.value,
			vepDb        = params.vepDb.value,
			filtervcf    = params.vepNonTCGAVcf.value,
			ref          = params.ref.value,
			oncotator    = params.oncotator.value,
			oncotator_db = params.oncotator_db.value,
			bcftools     = params.bcftools.value,
			tumor        = 'auto',
			genome       = params.genome.value,
			nthread      = 1,
			params       = Box()
		)
	)

@procfactory
def _pVcf2Plink():
	"""
	@name:
		pVcf2Plink
	@description:
		Convert vcf to plink binary files (.bed/.bim/.fam)
	@input:
		`infile:file`: The input vcf file, needs to be tabix indexed.
	@output:
		`outdir:dir`: The output directory containing the plink binary files
	@args:
		`plink` : The path to plink
		`params`: Command arguments for `plink`. Some pre-settings:
			- `vcf-half-call`      : `m`
			- `double-id`          : `True`
			- `vcf-filter`         : `True`
			- `vcf-idspace-to`     : `_`
			- `set-missing-var-ids`: `@_#`    # make sure no duplicate vars
				- if `$1`, `$2` ... included, this will run a extra process to set the var ids first
				- Since plink 1.x doesn't specify `$1` as ref, but the first one of all alleles in ASCII-sort order
				- Here `$1` will be bound to reference allele
			- `biallelic-only`     : `strict`
	@requires:
		`python:pyvcf`: to assign variant names (see `args.set-missing-var-ids`)
	"""
	pVcf2Plink             = Proc(desc = 'Convert vcf to plink binary files (.bed/.bim/.fam)')
	pVcf2Plink.input       = 'infile:file'
	pVcf2Plink.output      = 'outdir:dir:{{i.infile | fn2}}.plink'
	pVcf2Plink.args.plink  = params.plink.value
	pVcf2Plink.args.tabix  = params.tabix.value
	pVcf2Plink.args.params = Box({
		'vcf-half-call'      : 'm',
		'double-id'          : True,
		'vcf-filter'         : True,
		'vcf-idspace-to'     : '_',
		'set-missing-var-ids': '@_#', # may generate duplicate vars!
		'biallelic-only'     : 'strict'
	})
	pVcf2Plink.lang   = params.python.value
	pVcf2Plink.script = "file:scripts/vcf/pVcf2Plink.py"
	return pVcf2Plink

@procfactory
def _pVcfLiftover():
	"""
	@input:
		infile: The input vcf file
	@output:
		outfile: The output vcf file
		umfile:  The unmapped records
	@args:
		tool    (str) : Which tool to use
		picard  (str) : The path to picard
		lochain (file): The liftover chain file
		ref     (file): The reference genome
		mem     (str) : The memory to use
		tmpdir  (dir) : The temporary directory
		params  (Box) : The extra params for the tool
		bcftools(str) : Path to bcftools
			- Used to correct sample orders. `picard LiftoverVcf` sometimes swap samples.
	"""
	return Box(
		desc = 'Liftover VCF files',
		input = 'infile:file',
		output = 'outfile:file:{{i.infile | stem | stem}}.vcf, umfile:file:{{i.infile | stem | stem}}.unmapped.vcf',
		lang = params.python.value,
		args = Box(
			tool     = 'picard',
			bcftools = params.bcftools.value,
			picard   = params.picard.value,
			lochain  = params.lochain.value,
			ref      = params.ref.value,
			mem      = params.mem8G.value,
			tmpdir   = params.tmpdir.value,
			params   = Box()
		)
	)

@procfactory
def _pVcfAddChr():
	"""
	@name:
		pVcfAddChr
	@description:
		Add `chr` to records and contigs of vcf files.
	@args:
		`chr`: The prefix to add to each record.
	"""
	pVcfAddChr          = Proc(desc = 'Add `chr` to records and contigs of vcf files.')
	pVcfAddChr.input    = 'infile:file'
	pVcfAddChr.output   = 'outfile:file:{{i.infile | fn2}}.vcf'
	pVcfAddChr.args.chr = 'chr'
	pVcfAddChr.lang     = params.python.value
	pVcfAddChr.script   = "file:scripts/vcf/pVcfAddChr.py"
	return pVcfAddChr

@procfactory
def _pVcfCleanup():
	"""
	@input:
		infile: The input vcf file
	@output:
		outfile: The output vcf file
	@args:
		ref (file): The reference file
			- Require fai/dict file with it.
	"""
	return Box(
		desc   = 'Remove configs from vcf file according to the given reference',
		input  = 'infile:file',
		output = 'outfile:file:{{i.infile | stem | stem}}.vcf',
		lang   = params.python.value,
		args   = Box(ref = params.ref.value)
	)

@procfactory
def _pVcf2GTVcf():
	"""
	@name:
		pVcf2GTVcf
	@description:
		Keep only GT information for each sample.
	@input:
		`infile:file`: The input file
	@output:
		`outfile:file`: The output file, Default: `{{i.infile | fn2}}.vcf{{".gz" if args.gz else ""}}`
	@args:
		`tool`    : The tool to use, Default: `bcftools`
		`gz`      : Gzip the output file or not, Default: `False`
		`bcftools`: Path to bcftools, Default: `<params.bcftools>`
	"""
	pVcf2GTVcf               = Proc(desc = 'Keep only GT information for each sample.')
	pVcf2GTVcf.input         = 'infile:file'
	pVcf2GTVcf.output        = 'outfile:file:{{i.infile | fn2}}.vcf{{".gz" if args.gz else ""}}'
	pVcf2GTVcf.args.tool     = 'bcftools'
	pVcf2GTVcf.args.gz       = False
	pVcf2GTVcf.args.bcftools = params.bcftools.value
	pVcf2GTVcf.lang          = params.python.value
	pVcf2GTVcf.script        = "file:scripts/vcf/pVcf2GTVcf.py"
	return pVcf2GTVcf

@procfactory
def _pVcfSampleFilter():
	"""
	@name:
		pVcfSampleFilter
	@description:
		Keep or remove some samples from VCF file.
	@input:
		`infile:file` : The input file
		`samfile:file`: The file with sample names, one per line. Could be ignored, see `args.samples`
	@output:
		`outfile:file`: The output file, Default: `{{i.infile | fn2}}.vcf`
	@args:
		`bcftools`: Path to bcftools, Default: `<params.bcftools>`
		`samples`: Samples to filter, could be one of the followings, Default: `None`
			- A file with sample names, one per line.
			- A list with sample names
			- A string with sample names, separated by comma(,)
			- A string of lambda function to tell to keep current sample or not.
				- If this returns `True`, samples are added to the list, otherwise excluded.
				- Note that when `args.keep == False`, `True` samples will be removed.
			- `None`: use sample names from `i.samfile`
		`keep`: Keep the samples provided or remove them. Default: `True`
		`params`: Other parameters for `bcftools view`. Default: `Box(U = True)`
			- `U = True`: Exclude uncalled sites (genotypes of all samples are missing).
	"""
	pVcfSampleFilter               = Proc(desc = 'Keep or remove some samples from VCF file.')
	pVcfSampleFilter.input         = 'infile:file, samfile:file'
	pVcfSampleFilter.output        = 'outfile:file:{{i.infile | fn2}}.vcf'
	pVcfSampleFilter.args.keep     = True
	pVcfSampleFilter.args.samples  = None
	pVcfSampleFilter.args.params   = Box(U = True)
	pVcfSampleFilter.args.bcftools = params.bcftools.value
	pVcfSampleFilter.lang          = params.python.value
	pVcfSampleFilter.script        = "file:scripts/vcf/pVcfSampleFilter.py"
	return pVcfSampleFilter

@procfactory
def _pVcfSampleReplace():
	"""
	@name:
		pVcfSampleReplace
	@description:
		Replace sample names in VCF file
	@input:
		`infile:file` : The input file
		`samfile:file`: The file with sample names, one per line. Could be ignored, see `args.samples`
	@output:
		`outfile:file`: The output file, Default: `{{i.infile | fn2}}.vcf`
	@args:
		`bcftools`: Path to bcftools, Default: `<params.bcftools>`
		`samples`: New samples, could be one of the followings, Default: `None`
			- A file with new sample names, one per line.
			- A list with new sample names
			- A string with new sample names, separated by comma(,)
			- A string of lambda function to modify current sample names.
			- `None`: use sample names from `i.samfile`
		`nthread`: # threads used by `bcftools`
	"""
	pVcfSampleReplace               = Proc(desc = 'Replace sample names in VCF file')
	pVcfSampleReplace.input         = 'infile:file, samfile:file'
	pVcfSampleReplace.output        = 'outfile:file:{{i.infile | fn2}}.vcf'
	pVcfSampleReplace.args.bcftools = params.bcftools.value
	pVcfSampleReplace.args.samples  = None
	pVcfSampleReplace.args.nthread  = 0
	pVcfSampleReplace.lang          = params.python.value
	pVcfSampleReplace.script        = "file:scripts/vcf/pVcfSampleReplace.py"
	return pVcfSampleReplace

@procfactory
def _pVcf2GTMat():
	"""
	@name:
		pVcf2GTMat
	@description:
		Convert Vcf file to genotype matrix.
		Rownames are in the format of '<chr>_<pos>_<rs>_<ref>_<alt>'
	@input:
		`infile:file`: The input vcf file
	@output:
		`outfile:file`: the output filename. Default: `{{i.infile | fn2}}.gtmat`
	@args:
		`novel`: The snp name used if not mapped to any rsid. Default: `NOVEL`
		`useid`: Use the id in vcf file if possible. Default: `True`
		`dbsnp`: The dbsnp vcf file used to get the rsid. If not provided, will use `novel`
		`na`   : The value to replace missing genotypes.
		`bialt`: bi-allelic snps only. Default: `True`
	@requires:
		`pytabix`
		`pysam`
	"""
	pVcf2GTMat               = Proc(desc = 'Convert Vcf file to genotype matrix')
	pVcf2GTMat.input         = 'infile:file'
	pVcf2GTMat.output        = 'outfile:file:{{i.infile | fn2}}.gtmat.txt'
	pVcf2GTMat.args.novel    = 'NOVEL' # name. None to exclude variants without RSID
	pVcf2GTMat.args.useid    = True # use id in vcf file as possible
	pVcf2GTMat.args.dbsnp    = params.dbsnp_all.value
	pVcf2GTMat.args.tabix    = params.tabix.value
	pVcf2GTMat.args.samname  = None
	pVcf2GTMat.args.chrorder = params.chrorder.value
	pVcf2GTMat.args.bialt    = True # bi-allelic snps only
	pVcf2GTMat.args.mingt    = .5 # keep the variants with only the rate (mingt <= 1) or the number (mingt > 1)
	pVcf2GTMat.args.na       = 'NA'
	pVcf2GTMat.lang          = params.python.value
	pVcf2GTMat.script        = "file:scripts/vcf/pVcf2GTMat.py"
	return pVcf2GTMat

@procfactory
def _pVcfSort():
	"""
	@name:
		pVcfSort
	@description:
		Sort the vcf records
	@input:
		`infile:file`: The input file
	@output:
		`outfile:file`: The output file
	@args:
		`sortby`: Sort by what, Coordinates (coord) or names (name)? Default: `coord`
		`tool`  : The tool used to do the sort. Default: `sort` (linux command)
		`picard`: Path to picard.
		`tabix` : Path to tabix.
		`chrorder`: If sort by `args.sortby == 'coord'`, then records first sorted by `chrorder` then Coordinates.
	"""
	pVcfSort               = Proc(desc = 'Sort the vcf records')
	pVcfSort.input         = 'infile:file'
	pVcfSort.output        = 'outfile:file:{{i.infile | fn2}}.vcf'
	pVcfSort.args.sortby   = 'coord' # or name
	pVcfSort.args.tool     = 'sort' # picard
	pVcfSort.args.picard   = params.picard.value
	pVcfSort.args.tabix    = params.tabix.value
	pVcfSort.args.chrorder = params.chrorder.value
	pVcfSort.args.gsize    = params.gsize.value
	pVcfSort.args.nthread  = 1
	pVcfSort.lang          = params.python.value
	pVcfSort.script        = "file:scripts/vcf/pVcfSort.py"
	return pVcfSort

@procfactory
def _pVcfSubtract():
	"""
	@name:
		pVcfSubtract
	@description:
		Subtract one vcf file from another
	@input:
		`infile1:file`: The vcf file to be subtracted
		`infile2:file`: The background vcf file
	@output:
		`outfile:file`: The subtracted vcf file.
	@args:
		`header`  : Output header? Default: `True`
		`bychrom` : Split the vcf file by chromosomes, do subtraction and then merge them. Default: `False`
			- In case the vcf file is too big.
			- Requires both vcf files indexed (.tbi). If not they will be indexed there.
		`nthread` : # threads to use, only when `bychrom` is True. Default: `1`
		`tool`    : The tool to be used. Default: `mem` (or pyvcf/bedtools)
		`bedtools`: The path to bedtools.
		`tabix`   : The path to tabix.
		`any`     : Remove record in `infile1` with any overlap in `infile2`. Default: `True`
	"""
	pVcfSubtract               = Proc(desc = 'Subtract one vcf file from another')
	pVcfSubtract.input         = 'infile1:file, infile2:file'
	pVcfSubtract.output        = 'outfile:file:{{i.infile1 | fn2}}.subtracted.vcf'
	pVcfSubtract.args.bychrom  = False
	pVcfSubtract.args.nthread  = 1
	pVcfSubtract.args.header   = True
	pVcfSubtract.args.any      = True
	pVcfSubtract.args.tool     = 'mem'
	pVcfSubtract.args.tabix    = params.tabix.value
	pVcfSubtract.args.bedtools = params.bedtools.value
	pVcfSubtract.lang          = params.python.value
	pVcfSubtract.script        = "file:scripts/vcf/pVcfSubtract.py"
	return pVcfSubtract

@procfactory
def _pVcfExtract():
	"""
	@name:
		pVcfExtract
	@description:
		Extract variants from a VCF file by given regions
	@args:
		`tabix` : The path to tabix.
		`params`: Other parameters for `tabix`. Default: `Box(h = True, B = True)`
			- See `tabix --help`
	"""
	pVcfExtract              = Proc(desc = "Extract variants from a VCF file by given regions")
	pVcfExtract.input        = 'vcffile:file, regfile:file'
	pVcfExtract.output       = 'outfile:file:{{i.vcffile | fn2}}.extracted.vcf'
	pVcfExtract.args.tabix   = params.tabix.value
	pVcfExtract.args.params  = Box(h = True, B = True)
	pVcfExtract.lang         = params.python.value
	pVcfExtract.script       = "file:scripts/vcf/pVcfExtract.py"
	return pVcfExtract

@procfactory
def _pVcf2Pyclone():
	"""
	@name:
		pVcf2Pyclone
	"""
	pVcf2Pyclone        = Proc(desc = 'Generate PyClone input file for non-CN mutations')
	pVcf2Pyclone.input  = 'infile:file'
	pVcf2Pyclone.output = 'outfile:file:{{i.infile | bn}}.pyclone.txt'
	pVcf2Pyclone.lang   = params.python.value
	pVcf2Pyclone.script = "file:scripts/vcf/pVcf2Pyclone.py"
	return pVcf2Pyclone

@procfactory
def _pVcfFix():
	"""
	@input:
		infile: The input VCF file
	@output:
		outfile: The output fixed VCF file
	@args:
		ref: The reference genome.
			- fai/dict required to get valid contigs
		nthread (int): The number of threads used by openblas from numpy
		fixes: The issues to fix.
			- _inverse (bool): Inverse the fix switches. Only fix those ones with False.
				- The fixes with parameters (instead of True/False) will anyway be fixed.
			- clinvarLink (bool): Remove some clinvar links in INFO that are not well-formatted
			- addChr (bool): Try to add chr to chromosomes if not present.
			- addAF (bool): Try to add FORMAT/AF based on FORMAT/AD and FORMAT/DP.
			- tumorpos (bool|str|list): Try to put tumor sample before normal if it is Tumor-Normal paired VCF file. It could be:
				- False: to disable this fix
				- True: to match the file name to determine the tumor sample
				- str|list: A list (separated by comma in str) of tumor samples to match and identify the tumor sample
				- If it is not a 2-sample VCF file, this is disabled anyway.
			- headerInfo (bool|dict): Try to fix missing INFOs in header.
				- False: to disable this fix
				- True: to use `{ID: <info>, Number: 1, Type: String, Description: <info>.upper()}` to add INFO to header.
				- dict: to specify details to add INFO to header.
					- For example: `{COMMON: {Type: "Int", Description: "Whether it is a common SNP"}}`
			- headerFormat (bool|dict): Try to fix missing FORMATs in header.
				- Similar as headerInfo
			- headerFilter (bool|dict): Try to fix missing FILTERs in header.
				- False: to disable this fix
				- True: to use `<filter>.upper()` as description to add FILTER to header.
				- dict: to specify description for filters.
					- For example: `{LOWDEPTH: "Depth is lower than 50"}`
			- headerContig (bool): Add missing header contigs to the header.
				- `arg.ref` together with `<ref>.fai` or `<ref|stem>.dict` is recommended to get contig information
				- If `args.ref` is  not specified, a length of 999999999 will be used for missing contigs.
			- non_ref (bool): Fix the alternate alleles with `<NON_REF>`.
				- If `<NON_REF>` is the sole alternative allele, replace it with `.`
				- Otherwise remove it.
	"""
	return Box(
		desc   = 'Fix a bunch of format problems in vcf files',
		input  = 'infile:file',
		output = '''outfile:file:{{i.infile
			| ?.endswith(".gz") | : "_fixed.vcf.gz" | : "_fixed.vcf"
			| @prepend: stem(stem(i.infile))}}''',
		lang   = params.python.value,
		args   = Box(
			ref     = params.ref.value,
			nthread = 1,
			fixes   = Box(
				_inverse     = False,
				clinvarLink  = True,
				addChr       = True,
				addAF        = True,
				tumorpos     = True,
				headerInfo   = True,
				headerContig = True,
				headerFormat = True,
				headerFilter = True,
				non_ref      = True,
			)
		)
	)

@procfactory
def _pVcfFixGT():
	"""
	@name:
		pVcfFixGT
	@description:
		Fix genotypes in GVCF files.
		Some tools report genotype as 0/0 even there is no information for
		other FORMAT. This process replaces inappropriate 0/0 into ./.
	@input:
		`infile:file`: The input VCF file
	@output:
		`outfile:file`: The output VCF file, Default: `{{i.infile | bn}}`
	@args:
		`gz`: whether output (big) gzipped file. Default: `False`
	"""
	pVcfFixGT         = Proc(desc = 'Fix genotypes in GVCF files')
	pVcfFixGT.input   = 'infile:file'
	pVcfFixGT.output  = 'outfile:file:{{i.infile | bn}}'
	pVcfFixGT.args.gz = False
	pVcfFixGT.lang    = params.python.value
	pVcfFixGT.script  = "file:scripts/vcf/pVcfFixGT.py"
	return pVcfFixGT

@procfactory
def _pVcfStats():
	"""
	@input:
		infile: The input VCF file
		config: The configuration file for `vcfstats`
	@output:
		outdir: The output directory
	@args:
		vcfstats (str)      : Path to vcfstats.
		Rscript  (str)      : Path to Rscript.
		formula  (str/list) : Formulas to do the statistics.
		title    (str/list) : Title of each statistic.
		figtype  (str/list) : Type of figure for each statistic.
		passed   (bool)     : Whether Only take variants passed all filters in the statistics.
		region   (str/list) : Only take variants in region in the statistics, such as `chr1:1-1000`.
		regfile  (str)      : A bed file of regions.
		macro    (str)      : A macro for `vcfstats`.
		ggs      (str/list) : ggs expressions to modify each plot.
		devpars  (dict/list): Devpars for each plot.
	"""
	return Box(
		desc   = 'VCF statistics and plots using vcfstats',
		lang   = params.python.value,
		input  = 'infile:file, config:file',
		output = 'outdir:dir:{{i.infile | stem}}.vcfstats',
		args   = Box(
			vcfstats = params.vcfstats.value,
			Rscript  = params.Rscript.value,
			formula  = [],
			title    = [],
			figtype  = [],
			passed   = False,
			region   = [],
			regfile  = None,
			macro    = None,
			ggs      = [],
			devpars  = []))


