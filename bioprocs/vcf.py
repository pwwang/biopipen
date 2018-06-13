"""
A set of processes to generate/process vcf files
"""
from os import path
from glob import glob
from pyppl import Proc, Box
from . import params

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
	`filters`: The filters, should be a string of lambda function:
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
		- Remember if filters() returns True, record remained.
	`gz`     : Whether to gzip the output file. Default: False
	`keep`   : Whether to keep the filtered records. Default: True. (only for gatk, snpsift at filter step)
@requires:
	[`pyvcf`](https://github.com/jamescasbon/PyVCF)
"""
pVcfFilter              = Proc(desc = 'Filter records in vcf file.')
pVcfFilter.input        = "infile:file"
pVcfFilter.output       = "outfile:file:{{in.infile | fn}}.vcf{% if args.gz %}.gz{% endif %}"
pVcfFilter.args.fname   = 'pVcfFilter'
pVcfFilter.args.filters = None
pVcfFilter.args.gz      = False
pVcfFilter.args.keep    = True # only for gatk, snpsift at filter step
pVcfFilter.lang         = params.python.value
pVcfFilter.script       = "file:scripts/vcf/pVcfFilter.py"

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
pVcf.output         = "outfile:file:{{in.infile | fn}}.vcf{% if args.gz %}.gz{% endif %}"
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
	"outfile:file:{{in.infile | fn}}.{{args.tool}}/{{in.infile | fn}}.{{args.tool}}.vcf{% if args.gz %}.gz{% endif %}", 
	"outdir:dir:{{in.infile | fn}}.{{args.tool}}"
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
	`tool`:     The tool used to do extraction. Default: vcftools
	`vcftools`: The path of vcftools' vcf-subset
	`bcftools`: The path of bcftools, used to extract the sample names from input vcf file.
	`gatk`:     The path of gatk.
"""
pVcfSplit                     = Proc(desc = "Split multi-sample Vcf to single-sample Vcf files.")
pVcfSplit.input               = "infile:file, samples"
pVcfSplit.output              = "outdir:dir:{{in.infile | fn}}-individuals"
pVcfSplit.args.tool           = 'vcftools'
pVcfSplit.args.vcftools       = params.vcftools_subset.value
pVcfSplit.args.bcftools       = params.bcftools.value # used to extract samples
pVcfSplit.args.gatk           = params.gatk.value
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
pVcfMerge                     = Proc(desc = "Merge single-sample Vcf files to multi-sample Vcf file.")
pVcfMerge.input               = "infiles:files"
pVcfMerge.output              = "outfile:file:{{in.infiles[0] | fn}}_etc.vcf"

pVcfMerge.args.tool           = 'vcftools'
pVcfMerge.args.vcftools       = params.vcftools_merge.value
pVcfMerge.args.gatk           = params.gatk.value
pVcfMerge.args.ref            = params.ref.value # only for gatk
pVcfMerge.args.vep            = params.vep.value
pVcfMerge.args.nthread        = 1
pVcfMerge.lang                = params.python.value
pVcfMerge.script              = "file:scripts/vcf/pVcfMerge.py"

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
pVcf2Maf.output              = 'outfile:file:{{in.infile | fn | fn}}.maf'
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
	pVcf2GTMat
@description:
	Convert Vcf file to genotype matrix.
"""
pVcf2GTMat = Proc(desc = 'Convert Vcf file to genotype matrix')
pVcf2GTMat.input = 'infile:file'
pVcf2GTMat.output = 'outfile:file:{{in.infile | fn2}}.gtmat'
pVcf2GTMat.args.rnames = 'coord' # name
pVcf2GTMat.args.na = 'NA'
pVcf2GTMat.lang = params.python.value
pVcf2GTMat.script = "file:scripts/vcf/pVcf2GTMat.py"


