from pyppl import Proc, Box
from . import params, bashimport
from .utils import fs2name

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
	`bedtools`: The bedtools executable,                  default: "bedtools"
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
if [[ $(reffai_check {{args.ref | @append: ".fai" | squote}}) -eq 1 ]]; then
	reffai_do {{args.ref | squote}} {{args.ref | @append: ".fai" | squote}} {{args.samtools | squote}}
fi
'''
pBedGetfasta.lang   = params.python.value
pBedGetfasta.script = "file:scripts/bedtools/pBedGetfasta.py"


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
	`bedtools`: The bedtools executable, default: "bedtools"
	`params`:   Other parameters for `bedtools closest`, default: ""
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedClosest                 = Proc(desc = 'Find the closest elements')
pBedClosest.input           = "afile:file, bfile:file"
pBedClosest.output          = "outfile:file:{{i.afile | fn}}.closest.bt"
pBedClosest.args.bedtools   = params.bedtools.value
pBedClosest.args.params     = Box()
pBedClosest.envs.bashimport = bashimport
pBedClosest.script = '''
{{bashimport}} __init__.bash
{{args.bedtools}} closest $(cmdargs {{args.params | json | lambda x: __import__('json').loads(x) | quote}} '-' ' ') -a {{i.afile | squote}} -b {{i.bfile | squote}} > {{out.outfile | squote}}
'''

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
	`bedtools`: The bedtools executable, default: "bedtools"
	`params`:   Other parameters for `bedtools closest`, default: ""
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedClosest2                 = Proc(desc = 'Find the closest elements')
pBedClosest2.input           = "afile:file, bfiles:files"
pBedClosest2.output          = "outfile:file:{{i.afile | fn}}.closest.bt"
pBedClosest2.args.bedtools   = params.bedtools.value
pBedClosest2.args.params     = Box()
pBedClosest2.envs.bashimport = bashimport
pBedClosest2.script = '''
{{bashimport}} __init__.bash
{{args.bedtools}} closest $(cmdargs {{args.params | json | lambda x: __import__('json').loads(x) | quote}} '-' ' ') -a {{i.afile | squote}} -b {{i.bfiles | asquote}} > {{out.outfile | squote}}
'''

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
	`bin`:     The bedtools executable, default: "bedtools"
	`params`:  Other parameters for `bedtools flank`, default: ""
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
	`bedtools`: The bedtools executable, default: "bedtools"
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
	`bedtools`: The bedtools executable, default: "bedtools"
	`params`:   Other parameters for `bedtools intersect`, default: ""
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedIntersect2                     = Proc(desc = 'The wrapper of "bedtools intersect" with multiple input b files.')
pBedIntersect2.input               = "afile:file, bfiles:files"
pBedIntersect2.output              = "outfile:file:{{i.afile | fn}}.intersect2.bt"
pBedIntersect2.args.bedtools       = params.bedtools.value
pBedIntersect2.args.params         = Box()
pBedIntersect2.lang                = params.python.value
pBedIntersect2.script              = 'file:scripts/bedtools/pBedIntersect2.py'

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
	`bin`:     The bedtools executable, default: "bedtools"
	`informat`:The format of input file, whether is a "bed" file or "genome" size file. Default: "bed"
	`params`:  Other parameters for `bedtools makewindows`, default: ""
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedMakewindows = Proc(desc = 'Makes adjacent or sliding windows across a genome or BED file.')
pBedMakewindows.input         = "infile:file"
pBedMakewindows.output        = "outfile:file:{{i.infile | fn}}.window.bed"
pBedMakewindows.args.params   = Box()
pBedMakewindows.args.bedtools = params.bedtools.value
pBedMakewindows.args.intype   = 'bed'
pBedMakewindows.envs.bashimport = bashimport
pBedMakewindows.script = '''
{{bashimport}} __init__.bash
{{args.bedtools}} makewindows $(cmdargs {{args.params | json | lambda x: __import__('json').loads(x) | quote}} '-' ' ') {{args.intype | lambda x: '-b' if x=='bed' else '-g'}} {{i.infile | squote}} > {{out.outfile | squote}}
'''

# region pBedMerge2
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
	`bedtools`: The bedtools executable,               default: "bedtools"
	`params`  : Other parameters for `bedtools merge`, default: {}
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedMerge                     = Proc(desc = 'Merge regions in a bed file using `bedtools merge`.')
pBedMerge.input               = "infile:file"
pBedMerge.output              = "outfile:file:{{i.infile | fn}}.merged.bed"
pBedMerge.args.bedtools       = params.bedtools.value
pBedMerge.args.params         = Box()
pBedMerge.lang                = params.python.value
pBedMerge.script              = "file:scripts/bedtools/pBedMerge.py"
# endregion

# region pBedMerge2
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
	`bedtools`: The bedtools executable,               default: "bedtools"
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
	`bin`:     The bedtools executable, default: "bedtools"
	`params`:  Other parameters for `bedtools multiinter`, default: ""
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedMultiinter                 = Proc(desc = 'Identifies common intervals among multiple BED/GFF/VCF files.')
pBedMultiinter.input           = "infiles:files"
pBedMultiinter.output          = "outfile:file:{{i.infiles | fs2name}}.multiinter.bt"
pBedMultiinter.args.bedtools   = params.bedtools.value
pBedMultiinter.args.params     = Box()
pBedMultiinter.envs.fs2name    = fs2name
pBedMultiinter.envs.bashimport = bashimport
pBedMultiinter.script          = """
{{bashimport}} __init__.bash
{{args.bedtools}} multiinter $(cmdargs {{args.params | json | lambda x: __import__('json').loads(x) | quote}} '-' ' ') -i {{infiles | asquote}} > "{{outfile}}"
"""

# region pBedRandom
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
	`bedtools`: The bedtools executable,    default: "bedtools"
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

"""
@name:
	pBedShift
@description:
	`bedtools shift` will move each feature in a feature file by a user-defined number of bases. While something like this could be done with an awk '{OFS="\t" print $1,$2+<shift>,$3+<shift>}', bedtools shift will restrict the resizing to the size of the chromosome (i.e. no features before 0 or past the chromosome end).
@input:
	`infile:file`: The input file
	`gfile:file`:  The genome size file
@output:
	`outfile:file`: The result file
@args:
	`bin`:     The bedtools executable, default: "bedtools"
	`params`:  Other parameters for `bedtools shift`, default: ""
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedShift                 = Proc(desc = '`bedtools shift` will move each feature in a feature file by a user-defined number of bases.')
pBedShift.input           = "infile:file, gfile:file"
pBedShift.output          = "outfile:file:{{infile | fn}}.shifted.bed"
pBedShift.args.bedtools   = params.bedtools.value
pBedShift.args.params     = Box()
pBedShift.envs.bashimport = bashimport
pBedShift.script          = """
{{bashimport}} __init__.bash
{{args.bedtools}} shift $(cmdargs {{args.params | json | lambda x: __import__('json').loads(x) | quote}} '-' ' ') -i "{{infile}}" -g "{{gfile}}" > "{{outfile}}"
"""

"""
@name:
	pBedShuffle
@description:
	`bedtools shuffle` will randomly permute the genomic locations of a feature file among a genome defined in a genome file. One can also provide an exclusions BED/GFF/VCF file that lists regions where you do not want the permuted features to be placed. For example, one might want to prevent features from being placed in known genome gaps. shuffle is useful as a null basis against which to test the significance of associations of one feature with another.
@input:
	`infile:file`: The input file
	`gfile:file`:  The genome size file
@output:
	`outfile:file`: The result file
@args:
	`bin`:     The bedtools executable, default: "bedtools"
	`params`:  Other parameters for `bedtools shuffle`, default: ""
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedShuffle                 = Proc(desc = '`bedtools shuffle` will randomly permute the genomic locations of a feature file among a genome defined in a genome file.')
pBedShuffle.input           = "infile:file, gfile:file"
pBedShuffle.output          = "outfile:file:{{infile | fn}}.shuffled.bed"
pBedShuffle.args.bedtools   = params.bedtools.value
pBedShuffle.args.params     = Box()
pBedShuffle.envs.bashimport = bashimport
pBedShuffle.script          = """
{{bashimport}} __init__.bash
{{args.bedtools}} shuffle $(cmdargs {{args.params | json | lambda x: __import__('json').loads(x) | quote}} '-' ' ') -i "{{infile}}" -g "{{gfile}}" > "{{outfile}}"
"""

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
	`bin`:     The bedtools executable, default: "bedtools"
	`params`:  Other parameters for `bedtools subtract`, default: ""
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedSubtract                 = Proc(desc = '`bedtools subtract` searches for features in B that overlap A.')
pBedSubtract.input           = "afile:file, bfile:file"
pBedSubtract.output          = "outfile:file:{{afile | fn}}.subtracted.bed"
pBedSubtract.args.bedtools   = params.bedtools.value
pBedSubtract.args.params     = Box()
pBedSubtract.envs.bashimport = bashimport
pBedSubtract.script          = """
{{bashimport}} __init__.bash
{{args.bedtools}} subtract $(cmdargs {{args.params | json | lambda x: __import__('json').loads(x) | quote}} '-' ' ') -a "{{afile}}" -b "{{bfile}}" > {{outfile}}
"""

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
	`bin`:     The bedtools executable, default: "bedtools"
	`params`:  Other parameters for `bedtools window`, default: ""
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedWindow                 = Proc(desc = 'Similar to `bedtools intersect`, `window` searches for overlapping features in A and B.')
pBedWindow.input           = "afile:file, bfile:file"
pBedWindow.output          = "outfile:file:{{afile | fn}}.window.bed"
pBedWindow.args.bedtools   = params.bedtools.value
pBedWindow.args.params     = Box()
pBedWindow.envs.bashimport = bashimport
pBedWindow.script          = """
{{bashimport}} __init__.bash
{{args.bedtools}} window $(cmdargs {{args.params | json | lambda x: __import__('json').loads(x) | quote}} '-' ' ') -a "{{afile}}" -b "{{bfile}}" > "{{outfile}}"
"""

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
	`bedtools`: The bedtools executable, default: "bedtools"
	`params`:   Other parameters for `bedtools genomecov`, default: `Box(bg = True)`
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedGenomecov                 = Proc(desc = '`bedtools genomecov` computes histograms (default), per-base reports (-d) and BEDGRAPH (-bg) summaries of feature coverage (e.g., aligned sequences) for a given genome.')
pBedGenomecov.input           = "infile:file"
pBedGenomecov.output          = "outfile:file:{{infile | fn}}.genomecov.bt"
pBedGenomecov.args.bedtools   = params.bedtools.value
pBedGenomecov.args.params     = Box(bg = True)
pBedGenomecov.envs.bashimport = bashimport
pBedGenomecov.script          = """
{{bashimport}} __init__.bash
{{args.bedtools}} genomecov $(cmdargs {{args.params | json | lambda x: __import__('json').loads(x) | quote}} '-' ' ') -ibam "{{infile}}" > "{{outfile}}"
"""

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
	`bedtools`: The bedtools executable, default: "bedtools"
	`params`:   Other parameters for `bedtools cluster`, default: `Box()`
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedCluster = Proc(desc = 'Similar to merge, cluster report each set of overlapping or "book-ended" features in an interval file.')
pBedCluster.input           = "infile:file"
pBedCluster.output          = "outfile:file:{{infile | fn}}.clustered.bt"
pBedCluster.args.bedtools   = params.bedtools.value
pBedCluster.args.params     = Box()
pBedCluster.envs.bashimport = bashimport
pBedCluster.script          = """
{{bashimport}} __init__.bash
{{args.bedtools}} cluster $(cmdargs {{args.params | json | lambda x: __import__('json').loads(x) | quote}} '-' ' ') -i "{{infile}}" > "{{outfile}}"
"""

