"""
A set of processes to generate/process vcf files
"""

from pyppl import Proc, Box
from .utils import mem2, runcmd, buildref, checkref, helpers, plot
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
	`tool`: Which tool to use for filtering. Default: 'vcflib'
	`vcflib_vcffilter`: The path of vcffilter from vcflib. Default: 'vcffilter'
	`gatk`            : The path of gatk. Default: 'gatk'
	`snpsift`         : The path of snpsift. Default: 'SnpSift'
	`bcftools`        : The path of bcftools. Default: 'bcftools'
	`samtools`        : The path of samtools. Default: 'samtools' (used by gatk to generate reference index)
	`picard`          : The path of picard. Default: 'picard' (used by picard to generate reference dict) 
	`params`          : Other params of `tool`. Default: ""
	`mem`             : The memory to be used. Default: "4G" (only for snpsift and gatk)
	`gz`              : Whether to gzip the output file. Default: False
	`keep`            : Whether to keep the filtered records. Default: True. (only for gatk, snpsift at filter step)
	`ref`             : The path of reference. Default: "" (for gatk)
	`tmpdir`          : The path of tmpdir. Default: <system tmpdir> (only used by gatk and snpsift)
	`nthread`         : The path of Default: 1	
	`selectors`:   Select records by:
	  - type (snp, indel), sample genotypes (0, 1, 2), min genotype quality, filter (PASS, .)
	  - for example:
	    ```
		{"type": "snp", "genotype": {0: '0/0'}, "qual": 30}
		to select snps and whose genotype is '0/0' in 1st sample with quality >= 30
		{"genotype": {0: ['1/1', '0|1']}, "filter": ["PASS"]}
		to select records with PASS and genotype in 1st sample is '1/1' or '0/1'
		```
	`filters`:     Filters depend on the tool you use on INFO filelds
	  - format: `{"name1": "expression1", ...}`
	  - If a string is specified, will convert to `{<tool name>: <expression>}`
	  - Remember it filters OUT the records when ANY of the expression is true
@requires:
	[`pyvcf`](https://github.com/jamescasbon/PyVCF)
	[`gatk`](https://software.broadinstitute.org/gatk)
	[`bcftools`](http://www.htslib.org/doc/bcftools-1.2.html)
	[`snpsift`](http://snpeff.sourceforge.net/SnpSift.version_4_0.html)
	[`samtools`](https://github.com/samtools/samtools) if `args.ref` is not indexed, or bamutil is used for bam index file generation.
	[`picard`](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
"""
pVcfFilter                        = Proc(desc = 'Filter records in vcf file.')
pVcfFilter.input                  = "infile:file"
pVcfFilter.output                 = "outfile:file:{{in.infile | fn}}.vcf{% if args.gz %}.gz{% endif %}"
pVcfFilter.args.tool              = 'vcflib' # vcflib,    gatk, snpsift
pVcfFilter.args.vcflib            = params.vcflib_vcffilter.value
pVcfFilter.args.gatk              = params.gatk.value
pVcfFilter.args.snpsift           = params.snpsift.value
pVcfFilter.args.bcftools          = params.bcftools.value
pVcfFilter.args.samtools          = params.samtools.value
pVcfFilter.args.picard            = params.picard.value
pVcfFilter.args.selectors         = Box()
pVcfFilter.args.filters           = Box()
pVcfFilter.args.params            = Box()
pVcfFilter.args.mem               = params.mem4G.value
pVcfFilter.args.gz                = False
pVcfFilter.args.keep              = True # only for gatk, snpsift at filter step
pVcfFilter.args.ref               = params.ref.value # gatk
pVcfFilter.args.tmpdir            = params.tmpdir.value
pVcfFilter.args.nthread           = 1
pVcfFilter.tplenvs.mem2           = mem2.py
pVcfFilter.tplenvs.runcmd         = runcmd.py
pVcfFilter.tplenvs.buildrefIndex  = buildref.index.py
pVcfFilter.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pVcfFilter.lang                   = params.python.value
pVcfFilter.script                 = "file:scripts/vcf/pVcfFilter.py"

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
pVcfAnno.args.snpeffStats       = False
pVcfAnno.args.params            = Box()
pVcfAnno.args.mem               = params.mem8G.value
pVcfAnno.tplenvs.runcmd         = runcmd.py
pVcfAnno.tplenvs.mem2           = mem2.py
pVcfAnno.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pVcfAnno.beforeCmd              = """
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
pVcfSplit                        = Proc(desc = "Split multi-sample Vcf to single-sample Vcf files.")
pVcfSplit.input                  = "infile:file, samples"
pVcfSplit.output                 = "outdir:dir:{{in.infile | fn}}-individuals"
pVcfSplit.args.tool              = 'vcftools'
pVcfSplit.args.vcftools          = params.vcftools_subset.value
pVcfSplit.args.bcftools          = params.bcftools.value # used to extract samples
pVcfSplit.args.gatk              = params.gatk.value
pVcfSplit.args.ref               = params.ref.value # only for gatk
pVcfSplit.args.nthread           = 1
pVcfSplit.tplenvs.runcmd         = runcmd.py
pVcfSplit.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pVcfSplit.lang                   = params.python.value
pVcfSplit.script                 = "file:scripts/vcf/pVcfSplit.py"


