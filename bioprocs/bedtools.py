"""Utilities from bedtools"""
from pyppl import Proc
from diot import Diot
from .utils import fs2name
from . import params, proc_factory

pBedGetfasta = proc_factory(
	desc   = '`bedtools getfasta` extracts sequences from a FASTA file for each of the intervals defined in a BED file.',
	lang   = params.python.value,
	config = Diot(annotate = """
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
"""))
pBedGetfasta.input             = "infile:file"
pBedGetfasta.output            = "outfile:file:{{i.infile | fn}}.fa"
pBedGetfasta.args.samtools     = params.samtools.value
pBedGetfasta.args.bedtools     = params.bedtools.value
pBedGetfasta.args.params       = Diot(name = True)
pBedGetfasta.args.ref          = params.ref.value
pBedGetfasta.config.runcmd_pre = '''
{{"reference.bash" | bashimport}}
export samtools={{args.samtools | squote}}
reference fasta {{args.ref | squote}}
'''


pBedClosest = proc_factory(
	desc   = 'Find the closest elements',
	lang   = params.python.value,
	config = Diot(annotate = """
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
"""))
pBedClosest.input         = "afile:file, bfile:file"
pBedClosest.output        = "outfile:file:{{i.afile | fn}}.closest.bt"
pBedClosest.args.bedtools = params.bedtools.value
pBedClosest.args.params   = Diot()

pBedClosest2 = proc_factory(
	desc   = 'Find the closest elements',
	lang   = params.python.value,
	config = Diot(annotate = """
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
"""))
pBedClosest2.input         = "afile:file, bfiles:files"
pBedClosest2.output        = "outfile:file:{{i.afile | fn}}.closest.bt"
pBedClosest2.args.bedtools = params.bedtools.value
pBedClosest2.args.params   = Diot()

pBedFlank = proc_factory(
	desc   = 'Create two new flanking intervals for each interval in a BED file.',
	lang   = params.python.value,
	config = Diot(annotate = """
	@description:
		`bedtools flank` will create two new flanking intervals for each interval in a BED file. Note that flank will restrict the created flanking intervals to the size of the chromosome (i.e. no start < 0 and no end > chromosome size).
	@input:
		infile: The input file
	@output:
		outfile: The result file
	@args:
		bedtools (path): The bedtools executable
		params   (Diot) : Other parameters for `bedtools flank`
		gfile    (path): The genome size file
		extend   (bool): Whether extend the flanking regions to include the element itself or not.
			- For example: `chr1:100-200` with `args.params.b = 10` will extend to `chr1:90-210`
			- But if `args.extend = False`, it will be `chr1:90-100` and `chr1:200-210`
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""),
	input  = "infile:file",
	output = "outfile:file:{{i.infile | fn}}.flank.bed",
	args   = Diot(
		extend   = False,
		gsize    = params.gsize.value,
		params   = Diot(),
		bedtools = params.bedtools.value
	)
)

pBedIntersect = proc_factory(
	desc   = 'The wrapper of "bedtools intersect" with single input b file.',
	lang   = params.python.value,
	config = Diot(annotate = """
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
	"""))
pBedIntersect.input         = "afile:file, bfile:file"
pBedIntersect.output        = "outfile:file:{{i.afile | fn}}.intersect.bt"
pBedIntersect.args.bedtools = params.bedtools.value
pBedIntersect.args.params   = Diot()

pBedIntersect2 = proc_factory(
	desc   = 'The wrapper of "bedtools intersect" with multiple input b files.',
	lang   = params.python.value,
	config = Diot(annotate = """
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
	"""))
pBedIntersect2.input         = "afile:file, bfiles:files"
pBedIntersect2.output        = "outfile:file:{{i.afile | fn}}.intersect2.bt"
pBedIntersect2.args.bedtools = params.bedtools.value
pBedIntersect2.args.params   = Diot()

pBedMakewindows = proc_factory(
	desc   = 'Makes adjacent or sliding windows across a genome or BED file.',
	lang   = params.python.value,
	config = Diot(annotate = """
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
		`params`:   Other parameters for `bedtools makewindows`, default: `Diot()`
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""))
pBedMakewindows.input         = "infile:file"
pBedMakewindows.output        = "outfile:file:{{i.infile | fn}}.window.bed"
pBedMakewindows.args.params   = Diot()
pBedMakewindows.args.bedtools = params.bedtools.value
pBedMakewindows.args.intype   = 'bed'

# region pBedMerge2
pBedMerge = proc_factory(
	desc   = 'Merge regions in a bed file using `bedtools merge`.',
	lang   = params.python.value,
	config = Diot(annotate = """
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
"""))
pBedMerge.input         = "infile:file"
pBedMerge.output        = "outfile:file:{{i.infile | fn}}.merged.bed"
pBedMerge.args.bedtools = params.bedtools.value
pBedMerge.args.params   = Diot()
# endregion

# region pBedMerge2
pBedMerge2 = proc_factory(
	desc   = 'A multi-input file model of `pBedMerge`: Merge multiple input files.',
	lang   = params.python.value,
	config = Diot(annotate = """
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
"""))
pBedMerge2.input         = "infiles:files"
pBedMerge2.output        = "outfile:file:{{i.infiles | fs2name}}.merged.bed"
pBedMerge2.args.bedtools = params.bedtools.value
pBedMerge2.args.params   = Diot()
pBedMerge2.envs.fs2name  = fs2name
# endregion

pBedMultiinter = proc_factory(
	desc   = 'Identifies common intervals among multiple BED/GFF/VCF files.',
	lang   = params.python.value,
	config = Diot(annotate = """
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
"""))
pBedMultiinter.input         = "infiles:files"
pBedMultiinter.output        = "outfile:file:{{i.infiles | fs2name}}.multiinter.bt"
pBedMultiinter.args.bedtools = params.bedtools.value
pBedMultiinter.args.params   = Diot()
pBedMultiinter.envs.fs2name  = fs2name

# region pBedRandom
pBedRandom = proc_factory(
	desc   = 'Generate a random set of intervals in BED format.',
	lang   = params.python.value,
	config = Diot(annotate = """
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
	"""))
pBedRandom.input         = "l, n"
pBedRandom.output        = "outfile:file:random.L{{i.l}}.N{{i.n}}.S{{args.seed}}.bed"
pBedRandom.args.seed     = None
pBedRandom.args.bedtools = params.bedtools.value
pBedRandom.args.gsize    = params.gsize.value
# endregion

pBedShift = proc_factory(
	desc   = '`bedtools shift` will move each feature in a feature file by a user-defined number of bases.',
	lang   = params.python.value,
	config = Diot(annotate = """
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
	"""))
pBedShift.input         = "infile:file"
pBedShift.output        = "outfile:file:{{infile | fn}}.shifted.bed"
pBedShift.args.bedtools = params.bedtools.value
pBedShift.args.gsize    = params.gsize.value
pBedShift.args.params   = Diot()

pBedShuffle = proc_factory(
	desc = 'Randomly permute the genomic locations of a BED file using `bedtools shuffle`',
	lang   = params.python.value,
	config = Diot(annotate = """
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
"""))
pBedShuffle.input         = "infile:file"
pBedShuffle.output        = "outfile:file:{{i.infile | fn}}.shuffled.bed"
pBedShuffle.args.bedtools = params.bedtools.value
pBedShuffle.args.n        = 0
pBedShuffle.args.gsize    = params.gsize.value
pBedShuffle.args.seed     = False
pBedShuffle.args.params   = Diot()

pBedSubtract = proc_factory(
	desc = '`bedtools subtract` searches for features in B that overlap A.',
	lang   = params.python.value,
	config = Diot(annotate = """
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
"""))
pBedSubtract.input         = "afile:file, bfile:file"
pBedSubtract.output        = "outfile:file:{{afile | fn}}.subtracted.bed"
pBedSubtract.args.bedtools = params.bedtools.value
pBedSubtract.args.params   = Diot()

pBedWindow = proc_factory(
	desc = 'Similar to `bedtools intersect`, `window` searches for overlapping features in A and B.',
	lang   = params.python.value,
	config = Diot(annotate = """
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
"""))
pBedWindow.input         = "afile:file, bfile:file"
pBedWindow.output        = "outfile:file:{{afile | fn}}.window.bed"
pBedWindow.args.bedtools = params.bedtools.value
pBedWindow.args.params   = Diot()

pBedGenomecov = proc_factory(
	desc = '`bedtools genomecov` computes histograms',
	lang   = params.python.value,
	config = Diot(annotate = """
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
		`params`:   Other parameters for `bedtools genomecov`, default: `Diot(bg = True)`
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""))
pBedGenomecov.input         = "infile:file"
pBedGenomecov.output        = "outfile:file:{{infile | fn}}.genomecov.bt"
pBedGenomecov.args.bedtools = params.bedtools.value
pBedGenomecov.args.params   = Diot(bg = True)

pBedCluster = proc_factory(
	desc   = 'Similar to merge, cluster report each set of overlapping or "book-ended" features in an interval file.',
	lang   = params.python.value,
	config = Diot(annotate = """
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
		`params`:   Other parameters for `bedtools cluster`, default: `Diot()`
	@requires:
		[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
	"""))
pBedCluster.input         = "infile:file"
pBedCluster.output        = "outfile:file:{{infile | fn}}.clustered.bt"
pBedCluster.args.bedtools = params.bedtools.value
pBedCluster.args.params   = Diot()
