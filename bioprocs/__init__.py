VERSION = "0.0.1a"

import json
from os import path
from tempfile import gettempdir
from sys import modules, stderr
from pyppl import params

DEFAULTS = {
	# path
	"tmpdir"          : gettempdir(),
	"tmpdir.desc"     : "The temporary directory.",
	"dbsnp"           : "",
	"dbsnp.desc"      : "The dbsnp in VCF format.",
	"kallistoIdx"     : "",
	"kallistoIdx.desc": "The kallisto index file.",
	"ref"             : "",
	"ref.desc"        : "The reference genome.",
	"refgene"         : "",
	"refgene.desc"    : "The reference genes in GTF format.",


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
	"bamstats"                 : "bamstats",
	"bamstats.desc"            : "The path of bamstats.",
	"bamutil"                  : "bam",
	"bamutil.desc"             : "The path of bamutil.",
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
	"fastqc"                   : "fastqc",
	"fastqc.desc"              : "The path of fastqc.",
	"gatk"                     : "gatk",
	"gatk.desc"                : "The path of gatk.",
	"htseq_count"              : "htseq-count",
	"htseq_count.desc"         : "The path of htseq-ount.",
	"kallisto"                 : "kallisto",
	"kallisto.desc"            : "The path of kallisto.",
	"multiqc"                  : "multiqc",
	"multiqc.desc"             : "The path of multiqc.",
	"ngm"                      : "ngm",
	"ngm.desc"                 : "The path of ngm.",
	"picard"                   : "picard",
	"picard.desc"              : "The path of picard.",
	"platypus"                 : "platypus",
	"platypus.desc"            : "The path of platypus.",
	"sambamba"                 : "sambamba",
	"sambamba.desc"            : "The path of sambamba.",
	"samtools"                 : "samtools",
	"samtools.desc"            : "The path of samtools.",
	"skewer"                   : "skewer",
	"skewer.desc"              : "The path of skewer.",
	"snvsniffer"               : "SNVSniffer",
	"snvsniffer.desc"          : "The path of snvsniffer.",
	"somaticsniper"            : "bam-somaticsniper",
	"somaticsniper.desc"       : "The path of bam-somaticsniper.",
	"star"                     : "star",
	"star.desc"                : "The path of star.",
	"strelka_germ"             : "configureStrelkaGermlineWorkflow.py",
	"strelka_germ.desc"        : "The path of configureStrelkaGermlineWorkflow.py.",
	"strelka_soma"             : "configureStrelkaSomaticWorkflow.py",
	"strelka_soma.desc"        : "The path of configureStrelkaSomaticWorkflow.py.",
	"trimmomatic"              : "trimmomatic",
	"trimmomatic.desc"         : "The path of trimmomatic.",
	"vardict"                  : "vardict",
	"vardict.desc"             : "The path of vardict.",
	"virmid"                   : "vardict",
	"virmid.desc"              : "The path of vardict.",
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
	path.join (path.expanduser('~'), ".bioProcs"),   # values overwritten
	path.join (path.expanduser('~'), ".bioProcs.json")
]
for cfgfile in cfgfiles:
	if not path.exists(cfgfile):
		continue
	params.loadFile (cfgfile)