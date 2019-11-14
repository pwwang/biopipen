__version__ = '0.1.2'

import inspect
from pathlib import Path
from tempfile import gettempdir
from sys import executable, modules
from pyppl import Proc
from pyparam import Params

# open to R (reticulate) to get the path of r util scripts
HERE    = Path(__file__).resolve().parent
DEFAULTS = {
	# constants
	"dbsnpver"        : "150",
	"dbsnpver.desc"   : "The dbsnp version used in cruzdb to query SNP information.",
	"genome"          : "hg19",
	"genome.desc"     : "Commonly used genome assembly.",
	"ncbikey"         : "",
	"ncbikey.desc"    : "The NCBI API Key for E-Utils.",
	"chrorder"        : ",".join('chr' + str(x) for x in list(range(1,23)) + ['X', 'Y', 'M', 'MT']),
	"chrorder.desc"   : "Preferrable chromsome orders, used for sorting.",

	# path
	"affysnps"           : "",
	"affysnps.desc"      : "Affymetrix SNP 6.0 Positions, used for THetA2 mostly.",
	"annovarDb"          : "",
	"annovarDb.desc"     : "The path of database for Annovar.",
	"cachedir"           : str(Path.home() / 'cache'),
	"cachedir.desc"      : "The directory to cache query data.",
	"consvdir"           : "",
	"consvdir.desc"      : "The directory containing conservation scores in bigWig files.\nUse ucsc-wig2bigwig the original files are wigFix files.",
	"cytoband"           : "",
	"cytoband.desc"      : "The cytoband file from ucsc database.",
	"dbsnp"              : "",
	"dbsnp.desc"         : "The dbsnp common variants in VCF format.",
	"dbsnp_all"          : "",
	"dbsnp_all.desc"     : "The dbsnp all variants in VCF format.",
	"gep70"              : "",
	"gep70.desc"         : "The GEP70 genes (col1: 51up, col2:19down).",
	"gsize"              : "",
	"gsize.desc"         : "The chromsome size file.",
	"highld"             : "",
	"highld.desc"        : "High LD regions",
	"ipidb"              : "",
	"ipidb.desc"         : "The IPI xref database file.",
	"kallistoIdx"        : "",
	"kallistoIdx.desc"   : "The kallisto index file.",
	"lcsig"              : "",
	"lcsig.desc"         : "Lung cancer gene signature (see https://www.ncbi.nlm.nih.gov/pubmed/19118056)",
	"liftover"           : "liftOver",
	"liftover.desc"      : "The liftover binary file.",
	"lochain"            : "",
	"lochain.desc"       : "The liftover chain file.",
	"mcr"                : "",
	"mcr.desc"           : "Matlab MCR path.",
	"mutsig_cvrg"        : "",
	"mutsig_cvrg.desc"   : "Coverage table for MutSig.",
	"mutsig_cvrt"        : "",
	"mutsig_cvrt.desc"   : "Covariates table for MutSig.",
	"mutsig_mutdict"     : "",
	"mutsig_mutdict.desc": "Mutation type dictionary file for MutSig.",
	"mutsig_chrdir"      : "",
	"mutsig_chrdir.desc" : "Directory containing reference sequence files for MutSig.",
	"oncotator_db"       : "",
	"oncotator_db.desc"  : "The data directory for oncotator.",
	"ref"                : "",
	"ref.desc"           : "The reference genome.",
	"refgene"            : "",
	"refgene.desc"       : "The reference genes in GTF format.",
	"refexon"            : "",
	"refexon.desc"       : "The exons of reference genes in GTF format.",
	"rsmap_gwsnp6"       : "",
	"rsmap_gwsnp6.desc"  : "GenomeWideSNP_6 RSID SNP probe mapping file.",
	"snpeffDb"           : "",
	"snpeffDb.desc"      : "The path of database for snpEff.",
	"superfreq_res"      : "",
	"superfreq_res.desc" : "superFreq resource directory.",
	"tflist"             : "",
	"tflist.desc"        : "The TF list file with motif as 1st column and TF name as 2nd column.",
	"tfmotifs"           : "",
	"tfmotifs.desc"      : "The TF motifs in MEME format.",
	"tmpdir"             : gettempdir(),
	"tmpdir.desc"        : "The temporary directory.",
	"vepNonTCGAVcf"      : "",
	"vepNonTCGAVcf.desc" : "The path of ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz for VEP.",
	"vepDb"              : "",
	"vepDb.desc"         : "The path of database for VEP.",


	# memoreis
	"mem4G"      : "4G",
	"mem4G.desc" : "4GB memory.",
	"mem8G"      : "8G",
	"mem8G.desc" : "8GB memory.",
	"mem16G"     : "16G",
	"mem16G.desc": "16GB memory.",
	"mem24G"     : "24G",
	"mem24G.desc": "24GB memory.",
	"mem32G"     : "32G",
	"mem32G.desc": "32GB memory.",

	# tools
	"allfit"                   : "All-FIT.py",
	"allfit.desc"              : "Path to All-FIT.py",
	"annovar"                  : "annotate_variation.pl",
	"annovar.desc"             : "The path of annovar's annotate_variation.pl.",
	"annovar_convert"          : "convert2annovar.pl",
	"annovar_convert.desc"     : "The path of annovar's convert2annovar.pl.",
	"awk"                      : "awk",
	"awk.desc"                 : "The awk program.",
	"arsample"                 : "sample",
	"arsample.desc"            : "The tool `sample` by Alex Reynolds.",
	"bam_readcount"            : "bam-readcount",
	"bam_readcount.desc"       : "Program to generate metrics at single nucleotide positions",
	"bamstats"                 : "bamstats",
	"bamstats.desc"            : "The path of bamstats.",
	"bamutil"                  : "bam",
	"bamutil.desc"             : "The path of bamutil.",
	"bcftools"                 : "bcftools",
	"bcftools.desc"            : "The path of bcftools.",
	"bedops_sort"              : "sort-bed",
	"bedops_sort.desc"         : "The path of sort-bed of bedops.",
	"bedtools"                 : "bedtools",
	"bedtools.desc"            : "The path of bedtools.",
	"biobambam_bamsort"        : "bamsort",
	"biobambam_bamsort.desc"   : "The path of bamsort of biobambam.",
	"biobambam_bamtofastq"     : "bamtofastq",
	"biobambam_bamtofastq.desc": "The path of bamtofastq of biobambam.",
	"bowtie2"                  : "bowtie2",
	"bowtie2.desc"             : "The path of bamtofastq of bowtie2.",
	"bwa"                      : "bwa",
	"bwa.desc"                 : "The path of bamtofastq of bwa.",
	"bwtool"                   : "bwtool",
	"bwtool.desc"              : "The path of bwtool.",
	"cnvkit"                   : "cnvkit.py",
	"cnvkit.desc"              : "The path of cnvkit.py.",
	"cnvnator"                 : "cnvnator",
	"cnvnator.desc"            : "The path of cnvnator.",
	"cnvnator2vcf"             : "cnvnator2VCF.pl",
	"cnvnator2vcf.desc"        : "The path of cnvnator2VCF.pl.",
	"curl"                     : "curl",
	"curl.desc"                : "The path of command line tool curl.",
	"cutadapt"                 : "cutadapt",
	"cutadapt.desc"            : "The path of cutadapt.",
	"dot"                      : "dot",
	"dot.desc"                 : "The path of dot.",
	"dtoxog"                   : "dtoxog",
	"dtoxog.desc"              : "The path of dtoxog, a wrapper of D-ToxoG.",
	"dwgsim"                   : "dwgsim",
	"dwgsim.desc"              : "The path of dwgsim.",
	"elprep"                   : "elprep",
	"elprep.desc"              : "The path of elprep.",
	"fastqc"                   : "fastqc",
	"fastqc.desc"              : "The path of fastqc.",
	"fimo"                     : "fimo",
	"fimo.desc"                : "The path of fimo from MEME suite.",
	"gatk"                     : "gatk",
	"gatk.desc"                : "The path of gatk.",
	"gdc_client"               : "gdc-client",
	"gdc_client.desc"          : "The GDC client for TCGA data download.",
	"gistic"                   : "gistic2",
	"gistic.desc"              : "The path of gistic.",
	"htseq_count"              : "htseq-count",
	"htseq_count.desc"         : "The path of htseq-count.",
	"lichee"                   : "lichee",
	"lichee.desc"              : "The path of lichee.",
	"jvarkit"                  : "jvarkit",
	"jvarkit.desc"             : "The path of jvarkit.",
	"kallisto"                 : "kallisto",
	"kallisto.desc"            : "The path of kallisto.",
	"maf2vcf"                  : "maf2vcf.pl",
	"maf2vcf.desc"             : "The path of maf2vcf.pl.",
	"multiqc"                  : "multiqc",
	"multiqc.desc"             : "The path of multiqc.",
	"mutsig"                   : "mutsig",
	"mutsig.desc"              : "The path of run_MutSigCV.sh",
	"netmhc"                   : "netMHC",
	"netmhc.desc"              : "The path of netMHC.",
	"netmhcpan"                : "netMHCpan",
	"netmhcpan.desc"           : "The path of netMHCpan.",
	"netmhciipan"              : "netMHCIIpan",
	"netmhciipan.desc"         : "The path of netMHCIIpan.",
	"netmhccons"               : "netMHCCons",
	"netmhccons.desc"          : "The path of netMHCCons.",
	"ngm"                      : "ngm",
	"ngm.desc"                 : "The path of ngm.",
	"oncotator"                : "oncotator",
	"oncotator.desc"           : "The path of oncotator.",
	"picard"                   : "picard",
	"picard.desc"              : "The path of picard.",
	"platypus"                 : "platypus",
	"platypus.desc"            : "The path of platypus.",
	"plink"                    : "plink",
	"plink.desc"               : "Executable of plink.",
	"pyclone"                  : "PyClone",
	"pyclone.desc"             : "The path of PyClone",
	"sambamba"                 : "sambamba",
	"sambamba.desc"            : "The path of sambamba.",
	"samtools"                 : "samtools",
	"samtools.desc"            : "The path of samtools.",
	"schism"                   : "runSchism",
	"schism.desc"              : "The path of runSchism.",
	"skewer"                   : "skewer",
	"skewer.desc"              : "The path of skewer.",
	"snpeff"                   : "snpEff",
	"snpeff.desc"              : "The path of snpEff.",
	"snpsift"                  : "SnpSift",
	"snpsift.desc"             : "The path of SnpSift.",
	"snvsniffer"               : "SNVSniffer",
	"snvsniffer.desc"          : "The path of SNVSniffer.",
	"somaticsniper"            : "bam-somaticsniper",
	"somaticsniper.desc"       : "The path of bam-somaticsniper.",
	"star"                     : "STAR",
	"star.desc"                : "The path of star.",
	"strelka_germ"             : "configureStrelkaGermlineWorkflow.py",
	"strelka_germ.desc"        : "The path of configureStrelkaGermlineWorkflow.py.",
	"strelka_soma"             : "configureStrelkaSomaticWorkflow.py",
	"strelka_soma.desc"        : "The path of configureStrelkaSomaticWorkflow.py.",
	"tabix"                    : "tabix",
	"tabix.desc"               : "The path of tabix.",
	"theta2"                   : "RunTHetA.py",
	"theta2.desc"              : "The path of THetA2 executable.",
	"tomtom"                   : "tomtom",
	"tomtom.desc"              : "tomtom from MEME suite.",
	"topiary"                  : "topiary",
	"topiary.desc"             : "The path to topiary.",
	"trimmomatic"              : "trimmomatic",
	"trimmomatic.desc"         : "The path of trimmomatic.",
	"vardict"                  : "vardict",
	"vardict.desc"             : "The path of vardict.",
	"vcf2maf"                  : "vcf2maf.pl",
	"vcf2maf.desc"             : "The path of vcf2maf.pl.",
	"vcfanno"                  : "vcfanno",
	"vcfanno.desc"             : "The path of vcfanno.",
	"vcflib_vcffilter"         : "vcffilter",
	"vcflib_vcffilter.desc"    : "The path of vcflib vcffilter.",
	"vcfstats"                 : "vcfstats",
	"vcfstats.desc"            : "The path of vcfstats.",
	"vcftools"                 : "vcftools",
	"vcftools.desc"            : "The path of vcftools.",
	"vcftools_merge"           : "vcf-merge",
	"vcftools_merge.desc"      : "The path of vcftools' vcf-merge.",
	"vcftools_subset"          : "vcf-subset",
	"vcftools_subset.desc"     : "The path of vcftools' vcf-subset.",
	"vep"                      : "vep",
	"vep.desc"                 : "The path of vep.",
	"virmid"                   : "virmid",
	"virmid.desc"              : "The path of virmid.",
	"wgsim"                    : "wgsim",
	"wgsim.desc"               : "The path of wgsim.",
	"wandy"                    : "Wandy",
	"wandy.desc"               : "The path of Wandy.",

	# langs
	"python"      : "python",
	"python.desc" : "The path of python.",
	"Rscript"     : "Rscript",
	"Rscript.desc": "The path of Rscript.",
	"bash"        : "bash",
	"bash.desc"   : "The path of bash.",
}

params = Params()
params._load(DEFAULTS)
cfgfiles = [
	Path.home() / '.bioprocs.config', # values overwritten
	Path.home() / '.bioprocs.json',
]
for cfgfile in cfgfiles:
	if not cfgfile.exists():
		continue
	params._loadFile (cfgfile)

# lock the params in case the options are overwritten unintentionally.
params._locked = True

cachedir = Path(params.cachedir.value)
if not cachedir.exists():
	cachedir.mkdir()

rimport  = """
(function(...) {
	reticulate::use_python('%s', required = TRUE)
	bioprocs = reticulate::import('bioprocs')
	for (rfile in list(...)) {
		source(file.path(bioprocs$HERE, 'utils', rfile))
	}
})""" % executable

bashimport = """
function __bashimport__ () {
	for src in "$@"; do
		source "%s/utils/$src"
	done
}
__python__='%s'
__bashimport__""" % (HERE, executable)


EXT_MAP = {
	'Rscript': 'R',
	'python' : 'py',
	'python2': 'py',
	'python3': 'py',
}
FACTORY_CACHE = {}

def delefactory():
	"""The factory to give the delegator for modkit"""
	frame  = inspect.currentframe().f_back
	module = inspect.getmodule(frame)
	def delegator(proc):
		try:
			procfac =  module._mkenvs['_' + proc]
			if not procfac:
				raise KeyError
		except KeyError as exc:
			raise ImportError('No such process: {!r}'.format(proc)) from exc
		if not callable(procfac):
			raise ImportError('Wrong type of process factory: {!r} in module {!r}'.format(
				'_' + proc, module.__name__))
		return procfac()
	return delegator

def _procfactory(procfunc, pid, alias, mdname, doc):
	def factory():
		proc = procfunc()
		if isinstance(proc, dict):
			proc = Proc(**proc)
		proc.id = alias
		proc.props.origin = pid
		lang = Path(proc.lang).name
		ext  = '.' + EXT_MAP.get(lang, lang)
		if ext == '.R' and not proc.envs.get('rimport'):
			proc.envs.rimport = rimport
		if ext == '.bash' and not proc.envs.get('bashimport'):
			proc.envs.bashimport = bashimport
		script = HERE / 'scripts' / mdname / (pid + ext)
		if not proc.config.script and script.exists():
			proc.script = 'file:%s' % script
		report = HERE / 'reports' / mdname / (pid + '.md')
		if not proc.config.report and report.exists():
			proc.report = 'file:%s' % report
		return proc
	factory.__doc__ = doc
	return factory

def procfactory(procfunc):
	mdname = procfunc.__module__.split('.')[-1]
	pid    = procfunc.__name__.lstrip('_')
	module = modules[procfunc.__module__]
	args   = inspect.signature(procfunc).parameters
	alias  = args.get('alias')
	if alias:
		alias = alias.default
		if alias[0] == '_':
			alias = alias[1:]
		module._mkenvs['_' + alias] = _procfactory(procfunc, pid, alias, mdname, procfunc.__doc__)
	return _procfactory(procfunc, pid, pid, mdname, procfunc.__doc__)
