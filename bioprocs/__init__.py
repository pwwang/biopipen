VERSION = "0.0.1a1"

import json
from os import path, makedirs
from tempfile import gettempdir
from sys import modules, stderr, executable
from pyppl import params

# open to R (reticulate) to get the path of r util scripts
UTILS    = path.join(path.realpath(path.dirname(__file__)), 'utils')
DEFAULTS = {
	# constants
	"dbsnpver"        : "150",
	"dbsnpver.desc"   : "The dbsnp version used in cruzdb to query SNP information.",
	"genome"          : "hg19",
	"genome.desc"     : "Commonly used genome assembly.",
	"ncbikey"         : "",
	"ncbikey.desc"    : "The NCBI API Key for E-Utils.",

	# path
	"affysnps"           : "",
	"affysnps.desc"      : "Affymetrix SNP 6.0 Positions, used for THetA2 mostly.",
	"annovarDb"          : "",
	"annovarDb.desc"     : "The path of database for Annovar.",
	"cachedir"           : path.expanduser("~/.bioprocs/cache"),
	"cachedir.desc"      : "The directory to cache query data.",
	"consvdir"           : "",
	"consvdir.desc"      : "The directory containing conservation scores in bigWig files.\nUse ucsc-wig2bigwig the original files are wigFix files.",
	"cytoband"           : "",
	"cytoband.desc"      : "The cytoband file from ucsc database.",
	"dbsnp"              : "",
	"dbsnp.desc"         : "The dbsnp common variants in VCF format.",
	"dbsnp_all"          : "",
	"dbsnp_all.desc"     : "The dbsnp all variants in VCF format.",
	"exbaits"            : "",
	"exbaits.desc"       : "The exome sequencing bait file",
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
	"annovar"                  : "annotate_variation.pl",
	"annovar.desc"             : "The path of annovar's annotate_variation.pl.",
	"annovar_convert"          : "convert2annovar.pl",
	"annovar_convert.desc"     : "The path of annovar's convert2annovar.pl.",
	"awk"                      : "awk",
	"awk.desc"                 : "The awk program.",
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
	"htseq_count.desc"         : "The path of htseq-ount.",
	"kallisto"                 : "kallisto",
	"kallisto.desc"            : "The path of kallisto.",
	"multiqc"                  : "multiqc",
	"multiqc.desc"             : "The path of multiqc.",
	"mutsig"                   : "mutsig",
	"mutsig.desc"              : "The path of run_MutSigCV.sh",
	"ngm"                      : "ngm",
	"ngm.desc"                 : "The path of ngm.",
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
	"trimmomatic"              : "trimmomatic",
	"trimmomatic.desc"         : "The path of trimmomatic.",
	"vardict"                  : "vardict",
	"vardict.desc"             : "The path of vardict.",
	"vcf2maf"                  : "vcf2maf.pl",
	"vcf2maf.desc"             : "The path of vcf2maf.pl.",
	"vcflib_vcffilter"         : "vcffilter",
	"vcflib_vcffilter.desc"    : "The path of vcflib vcffilter.",
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

params.loadDict(DEFAULTS)
cfgfiles = [
	path.join (path.expanduser('~'), ".bioprocs.config"),   # values overwritten
	path.join (path.expanduser('~'), ".bioprocs.json")
]
for cfgfile in cfgfiles:
	if not path.exists(cfgfile):
		continue
	params.loadFile (cfgfile)

if not path.exists(params.cachedir.value):
	makedirs(params.cachedir.value)

rimport  = """
(function(...) {
	library(reticulate)
	bioprocs = import('bioprocs')
	for (rfile in list(...)) {
		source(file.path(bioprocs$UTILS, rfile))
	}
})"""

bashimport = """
function __bashimport__ () {
	for src in "$@"; do
		source "%s/$src"
	done
}
__python__='%s'
__bashimport__""" % (UTILS, executable)
