"""Utilities from bedtools"""
from pyppl import Proc, Box
from . import params, bashimport
from .utils import fs2name
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory()).exports('_p*')

@procfactory
def _pBedGetfasta():
	"""
	@name:
		pBedGetfasta
	@description:
		`bedtools getfasta` extracts sequences from a FASTA file for each of the intervals defined in a BED file.
	@input:
		`infile:file`: The input bed file
	@output:
		`outfile:file`: The generated fasta file
	@args:
		`ref`     : The fasta file
		`bedtools`: The bedtools executable,                  default: `<params.bedtools>`
		`params`  : Other parameters for `bedtools getfasta`, default: ""
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedGetfasta                 = Proc(desc = '`bedtools getfasta` extracts sequences from a FASTA file for each of the intervals defined in a BED file.')
	pBedGetfasta.input           = "infile:file"
	pBedGetfasta.output          = "outfile:file:{{i.infile | fn}}.fa"
	pBedGetfasta.args.samtools   = params.samtools.value
	pBedGetfasta.args.bedtools   = params.bedtools.value
	pBedGetfasta.args.params     = Box(name = True)
	pBedGetfasta.args.ref        = params.ref.value
	pBedGetfasta.envs.bashimport = bashimport
	pBedGetfasta.beforeCmd       = '''
	{{bashimport}} reference.bash
	export samtools={{args.samtools | squote}}
	reference fasta {{args.ref | squote}}
	'''
	pBedGetfasta.lang   = params.python.value
	pBedGetfasta.script = "file:scripts/bedtools/pBedGetfasta.py"
	return pBedGetfasta


@procfactory
def _pBedClosest():
	"""
	@name:
		pBedClosest
	@description:
		Similar to intersect, closest searches for overlapping features in A and B. In the event that no feature in B overlaps the current feature in A, closest will report the nearest (that is, least genomic distance from the start or end of A) feature in B. For example, one might want to find which is the closest gene to a significant GWAS polymorphism. Note that closest will report an overlapping feature as the closest that is, it does not restrict to closest non-overlapping feature. The following iconic cheatsheet summarizes the funcitonality available through the various optyions provided by the closest tool.
	@input:
		`afile:file`: The -a file
		`bfile:file`: The -b file
	@output:
		`outfile:file`: The result file
	@args:
		`bedtools`: The bedtools executable, default: `<params.bedtools>`
		`params`:   Other parameters for `bedtools closest`, default: ""
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedClosest               = Proc(desc = 'Find the closest elements')
	pBedClosest.input         = "afile:file, bfile:file"
	pBedClosest.output        = "outfile:file:{{i.afile | fn}}.closest.bt"
	pBedClosest.args.bedtools = params.bedtools.value
	pBedClosest.args.params   = Box()
	pBedClosest.script        = "file:scripts/bedtools/pBedClosest.py"
	return pBedClosest

@procfactory
def _pBedClosest2():
	"""
	@name:
		pBedClosest2
	@description:
		Multiple b-file version of pBedClosest
	@input:
		`afile:file`:   The -a file
		`bfiles:files`: The -b files
	@output:
		`outfile:file`: The result file
	@args:
		`bedtools`: The bedtools executable, default: `<params.bedtools>`
		`params`:   Other parameters for `bedtools closest`, default: ""
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedClosest2               = Proc(desc = 'Find the closest elements')
	pBedClosest2.input         = "afile:file, bfiles:files"
	pBedClosest2.output        = "outfile:file:{{i.afile | fn}}.closest.bt"
	pBedClosest2.args.bedtools = params.bedtools.value
	pBedClosest2.args.params   = Box()
	pBedClosest2.script        = "file:scripts/bedtools/pBedClosest2.py"
	return pBedClosest2

@procfactory
def _pBedFlank():
	"""
	@name:
		pBedFlank
	@description:
		`bedtools flank` will create two new flanking intervals for each interval in a BED file. Note that flank will restrict the created flanking intervals to the size of the chromosome (i.e. no start < 0 and no end > chromosome size).
	@input:
		`infile:file`:  The input file
		`gfile:file`:   The genome size file
	@output:
		`outfile:file`: The result file
	@args:
		`bedtools`: The bedtools executable, default: `<params.bedtools>`
		`params`:   Other parameters for `bedtools flank`, default: ""
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedFlank               = Proc(desc = 'Create two new flanking intervals for each interval in a BED file.')
	pBedFlank.input         = "infile:file"
	pBedFlank.output        = "outfile:file:{{i.infile | fn}}.flank.bed"
	pBedFlank.args.extend   = False
	pBedFlank.args.gsize    = params.gsize.value
	pBedFlank.args.params   = Box()
	pBedFlank.args.bedtools = params.bedtools.value
	pBedFlank.lang          = params.python.value
	pBedFlank.script        = "file:scripts/bedtools/pBedFlank.py"
	return pBedFlank

@procfactory
def _pBedIntersect():
	"""
	@name:
		pBedIntersect
	@description:
		By far, the most common question asked of two sets of genomic features is whether or not any of the features in the two sets overlap with one another. This is known as feature intersection. bedtools intersect allows one to screen for overlaps between two sets of genomic features. Moreover, it allows one to have fine control as to how the intersections are reported. bedtools intersect works with both BED/GFF/VCF and BAM files as input.
	@input:
		`afile:file` : The a file
		`bfile:file`: The b file
	@output:
		`outfile:file`: The result file
	@args:
		`bedtools`: The bedtools executable, default: `<params.bedtools>`
		`params`:   Other parameters for `bedtools intersect`, default: ""
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedIntersect                     = Proc(desc = 'The wrapper of "bedtools intersect" with single input b file.')
	pBedIntersect.input               = "afile:file, bfile:file"
	pBedIntersect.output              = "outfile:file:{{i.afile | fn}}.intersect.bt"
	pBedIntersect.args.bedtools       = params.bedtools.value
	pBedIntersect.args.params         = Box()
	pBedIntersect.lang                = params.python.value
	pBedIntersect.script              = 'file:scripts/bedtools/pBedIntersect.py'
	return pBedIntersect

@procfactory
def _pBedIntersect2():
	"""
	@name:
		pBedIntersect2
	@description:
		Multiple b-file version of pBedIntersect
	@input:
		`afile:file` : The a file
		`bfiles:files`: The b files
	@output:
		`outfile:file`: The result file
	@args:
		`bedtools`: The bedtools executable, default: `<params.bedtools>`
		`params`:   Other parameters for `bedtools intersect`, default: ""
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedIntersect2               = Proc(desc = 'The wrapper of "bedtools intersect" with multiple input b files.')
	pBedIntersect2.input         = "afile:file, bfiles:files"
	pBedIntersect2.output        = "outfile:file:{{i.afile | fn}}.intersect2.bt"
	pBedIntersect2.args.bedtools = params.bedtools.value
	pBedIntersect2.args.params   = Box()
	pBedIntersect2.lang          = params.python.value
	pBedIntersect2.script        = 'file:scripts/bedtools/pBedIntersect2.py'
	return pBedIntersect2

@procfactory
def _pBedMakewindows():
	"""
	@name:
		pBedMakewindows
	@description:
		Makes adjacent or sliding windows across a genome or BED file.
	@input:
		`infile:file`: The input file
	@output:
		`outfile:file`: The result file
	@args:
		`bedtools`: The bedtools executable, default: `<params.bedtools>`
		`intype`:   The format of input file, whether is a "bed" file or "genome" size file. Default: `bed`
		`params`:   Other parameters for `bedtools makewindows`, default: `Box()`
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedMakewindows               = Proc(desc = 'Makes adjacent or sliding windows across a genome or BED file.')
	pBedMakewindows.input         = "infile:file"
	pBedMakewindows.output        = "outfile:file:{{i.infile | fn}}.window.bed"
	pBedMakewindows.args.params   = Box()
	pBedMakewindows.args.bedtools = params.bedtools.value
	pBedMakewindows.args.intype   = 'bed'
	pBedMakewindows.script        = "file:scripts/bedtools/pBedMakewindows.py"
	return pBedMakewindows

# region pBedMerge2
@procfactory
def _pBedMerge():
	"""
	@name:
		pBedMerge
	@description:
		`bedtools merge` combines overlapping or book-ended features in an interval file into a single feature which spans all of the combined features.
	@input:
		`infile:file`: The input file
	@output:
		`outfile:file`: The result file
	@args:
		`bedtools`: The bedtools executable,               default: `<params.bedtools>`
		`params`  : Other parameters for `bedtools merge`, default: {}
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedMerge               = Proc(desc = 'Merge regions in a bed file using `bedtools merge`.')
	pBedMerge.input         = "infile:file"
	pBedMerge.output        = "outfile:file:{{i.infile | fn}}.merged.bed"
	pBedMerge.args.bedtools = params.bedtools.value
	pBedMerge.args.params   = Box()
	pBedMerge.lang          = params.python.value
	pBedMerge.script        = "file:scripts/bedtools/pBedMerge.py"
	# endregion
	return pBedMerge

# region pBedMerge2
@procfactory
def _pBedMerge2():
	"""
	@name:
		pBedMerge2
	@description:
		A multi-input file model of pBedMerge: Merge multiple input files.
	@input:
		`infiles:files`: The input files
	@output:
		`outfile:file`: The result file
	@args:
		`bedtools`: The bedtools executable,               default: `<params.bedtools>`
		`params`  : Other parameters for `bedtools merge`, default: {}
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedMerge2               = Proc(desc = 'A multi-input file model of `pBedMerge`: Merge multiple input files.')
	pBedMerge2.input         = "infiles:files"
	pBedMerge2.output        = "outfile:file:{{i.infiles | fs2name}}.merged.bed"
	pBedMerge2.args.bedtools = params.bedtools.value
	pBedMerge2.args.params   = Box()
	pBedMerge2.envs.fs2name  = fs2name
	pBedMerge2.lang          = params.python.value
	pBedMerge2.script        = "file:scripts/bedtools/pBedMerge2.py"
	# endregion
	return pBedMerge2

@procfactory
def _pBedMultiinter():
	"""
	@name:
		pBedMultiinter
	@description:
		Identifies common intervals among multiple BED/GFF/VCF files.
	@input:
		`infiles:files`: The input files
	@output:
		`outfile:file`: The result file
	@args:
		`bedtools`: The bedtools executable, default: `<params.bedtools>`
		`params`:   Other parameters for `bedtools multiinter`, default: ""
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedMultiinter               = Proc(desc = 'Identifies common intervals among multiple BED/GFF/VCF files.')
	pBedMultiinter.input         = "infiles:files"
	pBedMultiinter.output        = "outfile:file:{{i.infiles | fs2name}}.multiinter.bt"
	pBedMultiinter.args.bedtools = params.bedtools.value
	pBedMultiinter.args.params   = Box()
	pBedMultiinter.envs.fs2name  = fs2name
	pBedMultiinter.script        = "file:scripts/bedtools/pBedMultiinter.py"
	return pBedMultiinter

# region pBedRandom
@procfactory
def _pBedRandom():
	"""
	@name:
		pBedRandom
	@description:
		`bedtools random` will generate a random set of intervals in BED6 format. One can specify both the number (-n) and the size (-l) of the intervals that should be generated.
	@input:
		`gfile:file`: The genome size file
	@output:
		`outfile:file`: The result file
	@args:
		`bedtools`: The bedtools executable,    default: `<params.bedtools>`
		`seed`    : The seed for randomization, default: None
		`gsize`   : The chromsize file.
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedRandom               = Proc(desc = 'Generate a random set of intervals in BED format.')
	pBedRandom.input         = "l, n"
	pBedRandom.output        = "outfile:file:random.L{{i.l}}.N{{i.n}}.S{{args.seed}}.bed"
	pBedRandom.args.seed     = None
	pBedRandom.args.bedtools = params.bedtools.value
	pBedRandom.args.gsize    = params.gsize.value
	pBedRandom.lang          = params.python.value
	pBedRandom.script        = "file:scripts/bedtools/pBedRandom.py"
	# endregion
	return pBedRandom

@procfactory
def _pBedShift():
	"""
	@name:
		pBedShift
	@description:
		`bedtools shift` will move each feature in a feature file by a user-defined number of bases. While something like this could be done with an awk '{OFS="\t" print $1,$2+<shift>,$3+<shift>}', bedtools shift will restrict the resizing to the size of the chromosome (i.e. no features before 0 or past the chromosome end).
	@input:
		`infile:file`: The input file
	@output:
		`outfile:file`: The result file
	@args:
		`gsize`   : The genome size file. Default: `<params.gsize>`
		`bedtools`: The bedtools executable. Default: `<params.bedtools>`
		`params`  : Other parameters for `bedtools shift`. Default: ``
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedShift               = Proc(desc = '`bedtools shift` will move each feature in a feature file by a user-defined number of bases.')
	pBedShift.input         = "infile:file"
	pBedShift.output        = "outfile:file:{{infile | fn}}.shifted.bed"
	pBedShift.args.bedtools = params.bedtools.value
	pBedShift.args.gsize    = params.gsize.value
	pBedShift.args.params   = Box()
	pBedShift.script        = "file:scripts/bedtools/pBedShift.py"
	return pBedShift

@procfactory
def _pBedShuffle():
	"""
	@name:
		pBedShuffle
	@description:
		`bedtools shuffle` will randomly permute the genomic locations of a feature file among a genome defined in a genome file. One can also provide an exclusions BED/GFF/VCF file that lists regions where you do not want the permuted features to be placed. For example, one might want to prevent features from being placed in known genome gaps. shuffle is useful as a null basis against which to test the significance of associations of one feature with another.
	@input:
		`infile:file`: The input file
	@output:
		`outfile:file`: The result file
	@args:
		`bedtools`: The bedtools executable, default: `<params.bedtools>`
		`params`  : Other parameters for `bedtools shuffle`, default: ""
		`gsize`   : The chromsize file. Default: `params.gsize`
		`n`       : Only return top `n` records (act like sampling). Default: `0`
			- `0`/`None` will return all records
			- `n<1` will return proportion of the records
			- `n>=1` will return top `n` records
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedShuffle               = Proc(desc = 'Randomly permute the genomic locations of a BED file using `bedtools shuffle`')
	pBedShuffle.input         = "infile:file"
	pBedShuffle.output        = "outfile:file:{{i.infile | fn}}.shuffled.bed"
	pBedShuffle.args.bedtools = params.bedtools.value
	pBedShuffle.args.n        = 0
	pBedShuffle.args.gsize    = params.gsize.value
	pBedShuffle.args.seed     = False
	pBedShuffle.args.params   = Box()
	pBedShuffle.lang          = params.python.value
	pBedShuffle.script        = "file:scripts/bedtools/pBedShuffle.py"
	return pBedShuffle

@procfactory
def _pBedSubtract():
	"""
	@name:
		pBedSubtract
	@description:
		`bedtools subtract` searches for features in B that overlap A. If an overlapping feature is found in B, the overlapping portion is removed from A and the remaining portion of A is reported. If a feature in B overlaps all of a feature in A, the A feature will not be reported.
	@input:
		`afile:file`: The a file
		`bfile:file`: The b file
	@output:
		`outfile:file`: The result file
	@args:
		`bedtools`: The bedtools executable, default: `<params.bedtools>`
		`params`:   Other parameters for `bedtools subtract`, default: ""
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedSubtract               = Proc(desc = '`bedtools subtract` searches for features in B that overlap A.')
	pBedSubtract.input         = "afile:file, bfile:file"
	pBedSubtract.output        = "outfile:file:{{afile | fn}}.subtracted.bed"
	pBedSubtract.args.bedtools = params.bedtools.value
	pBedSubtract.args.params   = Box()
	pBedSubtract.script        = "file:scripts/bedtools/pBedSubtract.py"
	return pBedSubtract

@procfactory
def _pBedWindow():
	"""
	@name:
		pBedWindow
	@description:
		Similar to `bedtools intersect`, `window` searches for overlapping features in A and B. However, window adds a specified number (1000, by default) of base pairs upstream and downstream of each feature in A. In effect, this allows features in B that are near features in A to be detected.
	@input:
		`afile:file`: The a file
		`bfile:file`: The b file
	@output:
		`outfile:file`: The result file
	@args:
		`bedtools`: The bedtools executable, default: `<params.bedtools>`
		`params`:   Other parameters for `bedtools window`, default: ""
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedWindow               = Proc(desc = 'Similar to `bedtools intersect`, `window` searches for overlapping features in A and B.')
	pBedWindow.input         = "afile:file, bfile:file"
	pBedWindow.output        = "outfile:file:{{afile | fn}}.window.bed"
	pBedWindow.args.bedtools = params.bedtools.value
	pBedWindow.args.params   = Box()
	pBedWindow.script        = "file:scripts/bedtools/pBedWindow.py"
	return pBedWindow

@procfactory
def _pBedGenomecov():
	"""
	@name:
		pBedGenomecov
	@description:
		`bedtools genomecov` computes histograms (default), per-base reports (-d) and BEDGRAPH (-bg) summaries of feature coverage (e.g., aligned sequences) for a given genome.
		NOTE: only bam file input implemented here.
	@input:
		`infile:file`: The bam file
	@output:
		`outfile:file`: The result file
	@args:
		`bedtools`: The bedtools executable, default: `<params.bedtools>`
		`params`:   Other parameters for `bedtools genomecov`, default: `Box(bg = True)`
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedGenomecov               = Proc(desc = '`bedtools genomecov` computes histograms (default), per-base reports (-d) and BEDGRAPH (-bg) summaries of feature coverage (e.g., aligned sequences) for a given genome.')
	pBedGenomecov.input         = "infile:file"
	pBedGenomecov.output        = "outfile:file:{{infile | fn}}.genomecov.bt"
	pBedGenomecov.args.bedtools = params.bedtools.value
	pBedGenomecov.args.params   = Box(bg = True)
	pBedGenomecov.script        = "file:scripts/bedtools/pBedGenomecov.py"
	return pBedGenomecov

@procfactory
def _pBedCluster():
	"""
	@name:
		pBedCluster
	@description:
		Similar to merge, cluster report each set of overlapping or "book-ended" features in an interval file. In contrast to merge, cluster does not flatten the cluster of intervals into a new meta-interval; instead, it assigns an unique cluster ID to each record in each cluster. This is useful for having fine control over how sets of overlapping intervals in a single interval file are combined.
	@input:
		`infile:file`: The input file
	@output:
		`outfile:file`: The output file with cluster id for each record
	@args:
		`bedtools`: The bedtools executable, default: `<params.bedtools>`
		`params`:   Other parameters for `bedtools cluster`, default: `Box()`
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""
	pBedCluster               = Proc(desc = 'Similar to merge, cluster report each set of overlapping or "book-ended" features in an interval file.')
	pBedCluster.input         = "infile:file"
	pBedCluster.output        = "outfile:file:{{infile | fn}}.clustered.bt"
	pBedCluster.args.bedtools = params.bedtools.value
	pBedCluster.args.params   = Box()
	pBedCluster.lang          = params.python.value
	pBedCluster.script        = "file:scripts/bedtools/pBedCluster.py"
	return pBedCluster

