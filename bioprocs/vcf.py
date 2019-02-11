# A set of processes to generate/process vcf files
from os import path
from glob import glob
from pyppl import Proc, Box
from . import params
from .utils import fs2name

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
	`readerops`: A lambda function (must be quoted) to manipulate the reader (vcf.Reader instance)
	`recordops`: A lambda function (must be quoted) to manipulate the record (vcf.Record instance)
	`gz`: Gzip the ouput file
"""
pVcf                = Proc(desc = 'Motify vcf file using pyvcf')
pVcf.input          = "infile:file"
pVcf.output         = "outfile:file:{{i.infile | fn}}.vcf{% if args.gz %}.gz{% endif %}"
pVcf.args.helper    = ''
pVcf.args.readerops = 'lambda x: None'
pVcf.args.recordops = 'lambda x, y = None: None'
pVcf.args.gz        = False
pVcf.lang           = params.python.value
pVcf.script         = "file:scripts/vcf/pVcf.py"

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
	`dbpath`:          The path of database for each tool. Required by 'annovar' and 'vep'
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
pVcfAnno.args.tool            = 'snpeff'
pVcfAnno.args.snpeff          = params.snpeff.value
pVcfAnno.args.vep             = params.vep.value
pVcfAnno.args.gz              = False
pVcfAnno.args.annovar         = params.annovar.value
pVcfAnno.args.annovar_convert = params.annovar_convert.value
pVcfAnno.args.genome          = params.genome.value
pVcfAnno.args.tmpdir          = params.tmpdir.value
pVcfAnno.args.dbpath          = Box({
	'snpeff' : params.snpeffDb.value,
	'annovar': params.annovarDb.value,
	'vep'    : params.vepDb.value
})
pVcfAnno.args.snpeffStats    = False
pVcfAnno.args.params         = Box()
pVcfAnno.args.mem            = params.mem8G.value
pVcfAnno.beforeCmd           = """
# check dbpath
dbpath=$({{proc.lang}} -c "print {{args.dbpath}}['{{args.tool}}']")
if [[ ! -e "$dbpath" ]]; then
	echo "You have to specify valid db path for tool: {{args.tool}}" 1>&2 
	echo "  - For vep: /path/to/cache" 1>&2
	echo "  - For snpEff: /path/to/datadir" 1>&2
	echo "  - For annovar: /path/to/db" 1>&2
	exit 1
fi
"""
pVcfAnno.lang                 = params.python.value
pVcfAnno.script               = "file:scripts/vcf/pVcfAnno.py"

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
	`tool`:     The tool used to do extraction. Default: vcftools (gatk, awk)
	`vcftools`: The path of vcftools' vcf-subset
	`bcftools`: The path of bcftools, used to extract the sample names from input vcf file.
	`gatk`:     The path of gatk.
"""
pVcfSplit                     = Proc(desc = "Split multi-sample Vcf to single-sample Vcf files.")
pVcfSplit.input               = "infile:file, samples"
pVcfSplit.output              = "outdir:dir:{{i.infile | fn}}-individuals"
pVcfSplit.args.tool           = 'vcftools'
pVcfSplit.args.vcftools       = params.vcftools_subset.value
pVcfSplit.args.bcftools       = params.bcftools.value # used to extract samples
pVcfSplit.args.gatk           = params.gatk.value
pVcfSplit.args.awk            = params.awk.value
pVcfSplit.args.ref            = params.ref.value # only for gatk
pVcfSplit.args.params         = Box()
pVcfSplit.args.nthread        = 1
pVcfSplit.lang                = params.python.value
pVcfSplit.script              = "file:scripts/vcf/pVcfSplit.py"

"""
@name:
	pVcfMerge
@description:
	Merge single-sample Vcf files to multi-sample Vcf file.
@input:
	`infiles:files`: The input vcf files
	`outfile:dir`:  The output multi-sample vcf.
@args:
	`tool`:     The tool used to do extraction. Default: vcftools
	`vcftools`: The path of vcftools' vcf-subset
	`bcftools`: The path of bcftools, used to extract the sample names from input vcf file.
	`gatk`:     The path of gatk.
"""
pVcfMerge               = Proc(desc = "Merge single-sample Vcf files to multi-sample Vcf file.")
pVcfMerge.input         = "infiles:files"
pVcfMerge.output        = "outfile:file:{{i.infiles | fs2name}}.vcf"
pVcfMerge.args.tool     = 'vcftools'
pVcfMerge.args.vcftools = params.vcftools_merge.value
pVcfMerge.args.gatk     = params.gatk.value
pVcfMerge.args.params   = Box()
pVcfMerge.args.tabix    = params.tabix.value
pVcfMerge.args.ref      = params.ref.value # only for gatk
pVcfMerge.args.nthread  = 1
pVcfMerge.envs.fs2name  = fs2name
pVcfMerge.lang          = params.python.value
pVcfMerge.script        = "file:scripts/vcf/pVcfMerge.py"

"""
@name:
	pVcf2Maf
@description:
	Convert Vcf file to Maf file
@input:
	`infile:file` : The input vcf file
		- see `args.somatic`
@output:
	`outfile:file`: The output maf file
@args:
	`tool`     : Which tool to use. Default: vcf2maf
	`vcf2maf`  : The path of vcf2maf.pl
	`vep`      : The path of vep
	`vepDb`    : The path of database for vep
	`filtervcf`: The filter vcf. Something like: ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
	`ref`      : The reference genome
	`nthread`  : Number of threads used to extract samples. Default: 1
	`tumor1st` : Whether tumor sample comes first. Default: `True`
	`bcftools` : Path to bcftools used to extract sample names.
	`vcftools` : Path to vcftools used to split vcf.
	`samfunc`  : A lambda function used to deduce sample names from file name.
	`somatic`  : Whether input vcf file is a somatic mutation file. Default: False
		- somatic mutation vcf file can only have one sample TUMOR, or two samples, TUMOR and NORMAL, but will be considered as single sample.
		- otherwise, multiple samples are supported in the input vcf file. Tumor id will be sample name for each sample, normal id will be NORMAL.
"""
pVcf2Maf                     = Proc(desc = 'Convert Vcf file to Maf file.')
pVcf2Maf.input               = 'infile:file'
pVcf2Maf.output              = 'outfile:file:{{i.infile | fn | fn}}.maf'
pVcf2Maf.args.tool           = 'vcf2maf'
pVcf2Maf.args.vcf2maf        = params.vcf2maf.value
pVcf2Maf.args.vep            = params.vep.value
pVcf2Maf.args.vepDb          = params.vepDb.value
pVcf2Maf.args.filtervcf      = params.vepNonTCGAVcf.value
pVcf2Maf.args.ref            = params.ref.value
pVcf2Maf.args.bcftools       = params.bcftools.value
pVcf2Maf.args.vcftools       = params.vcftools_subset.value
pVcf2Maf.args.samfunc        = None
pVcf2Maf.args.tumor1st       = True
pVcf2Maf.args.somatic        = False
pVcf2Maf.args.nthread        = 1
pVcf2Maf.args.params         = Box()
pVcf2Maf.lang                = params.python.value
pVcf2Maf.script              = "file:scripts/vcf/pVcf2Maf.py"

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
		- `biallelic-only`     : `strict`
"""
pVcf2Plink             = Proc(desc = 'Convert vcf to plink binary files (.bed/.bim/.fam)')
pVcf2Plink.input       = 'infile:file'
pVcf2Plink.output      = 'outdir:dir:{{i.infile | fn2}}.plink'
pVcf2Plink.args.plink  = params.plink.value
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

"""
@name:
	pVcfLiftover
@description:
	Lift over vcf files.
@input:
	`infile:file`: The input vcf file.
@output:
	`outfile:file`: The output vcf file.
	`umfile:file`:  The unmapped records
@args:
	`tool`:    Which tool to use. Default: `picard`
	`picard`:  The path to picard
	`lochain`: The liftover chain file
	`ref`:     The reference genome
	`mem`:     The memory to use
	`tmpdir`:  The temporary directory
	`params`:  The extra params.
"""
pVcfLiftover              = Proc(desc = 'Lift over vcf files.')
pVcfLiftover.input        = 'infile:file'
pVcfLiftover.output       = 'outfile:file:{{i.infile | fn}}.vcf, umfile:file:{{i.infile | fn}}.unmapped.vcf'
pVcfLiftover.args.tool    = 'picard'
pVcfLiftover.args.picard  = params.picard.value
pVcfLiftover.args.lochain = params.lochain.value
pVcfLiftover.args.ref     = params.ref.value
pVcfLiftover.args.mem     = params.mem8G.value
pVcfLiftover.args.tmpdir  = params.tmpdir.value
pVcfLiftover.args.params  = Box()
pVcfLiftover.lang         = params.python.value
pVcfLiftover.script       = "file:scripts/vcf/pVcfLiftover.py"

"""
@name:
	pVcfAddChr
@description:
	Add `chr` to records of vcf files.
@args:
	`chr`: The prefix to add to each record.
"""
pVcfAddChr          = Proc(desc = 'Add `chr` to records of vcf files.')
pVcfAddChr.input    = 'infile:file'
pVcfAddChr.output   = 'outfile:file:{{i.infile | fn}}.vcf'
pVcfAddChr.args.chr = 'chr'
pVcfAddChr.lang     = params.python.value
pVcfAddChr.script   = "file:scripts/vcf/pVcfAddChr.py"

"""
@name:
	pVcfCleanup
@description:
	Remove configs from vcf file according to the given reference.
@input:
	`infile:file`: The input vcf file
@output:
	`outfile:file`: The output vcf file. Default: `{{i.infile | fn}}.vcf`
@args:
	`ref`: The reference file
"""
pVcfCleanup          = Proc(desc = 'Remove configs from vcf file according to the given reference.')
pVcfCleanup.input    = 'infile:file'
pVcfCleanup.output   = 'outfile:file:{{i.infile | fn}}.vcf'
pVcfCleanup.args.ref = params.ref.value # required dict or fai file with it
pVcfCleanup.lang     = params.python.value
pVcfCleanup.script   = "file:scripts/vcf/pVcfCleanup.py"

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
	`useid`: Use the id in vcf file is possible. Default: `True`
	`dbsnp`: The dbsnp vcf file used to get the rsid. If not provided, will use `novel`
	`na`   : The value to replace missing genotypes.
	`bialt`: bi-allelic snps only. Default: `True`
@requires:
	`pytabix`
	`pysam`
"""
pVcf2GTMat            = Proc(desc = 'Convert Vcf file to genotype matrix')
pVcf2GTMat.input      = 'infile:file'
pVcf2GTMat.output     = 'outfile:file:{{i.infile | fn2}}.gtmat'
pVcf2GTMat.args.novel = 'NOVEL' # name
pVcf2GTMat.args.useid = True # use id in vcf file as possible
pVcf2GTMat.args.dbsnp = params.dbsnp_all.value
pVcf2GTMat.args.bialt = True # bi-allelic snps only
pVcf2GTMat.args.na    = 'NA'
pVcf2GTMat.lang       = params.python.value
pVcf2GTMat.script     = "file:scripts/vcf/pVcf2GTMat.py"

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
	`header`: Output header? Default: `True`
	`by`    : Sort by what, Coordinates (coord) or names (name)? Default: `coord`
	`tool`  : The tool used to do the sort. Default: `sort` (linux command)
"""
pVcfSort               = Proc(desc = 'Sort the vcf records')
pVcfSort.input         = 'infile:file'
pVcfSort.output        = 'outfile:file:{{i.infile | fn2}}.vcf'
pVcfSort.args.header   = True
pVcfSort.args.by       = 'coord' # or name
pVcfSort.args.tool     = 'sort' # or bedtools
pVcfSort.lang          = params.python.value
pVcfSort.script        = "file:scripts/vcf/pVcfSort.py"

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


pVcf2Pyclone        = Proc(desc = 'Generate PyClone input file for non-CN mutations')
pVcf2Pyclone.input  = 'infile:file'
pVcf2Pyclone.output = 'outfile:file:{{i.infile | bn}}.pyclone.txt'
pVcf2Pyclone.lang   = params.python.value
pVcf2Pyclone.script = "file:scripts/vcf/pVcf2Pyclone.py"

