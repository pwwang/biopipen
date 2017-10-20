from pyppl import Proc, Box
from . import params
from .utils import runcmd, helpers, buildref, checkref

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
pBedGetfasta                        = Proc(desc = '`bedtools getfasta` extracts sequences from a FASTA file for each of the intervals defined in a BED file.')
pBedGetfasta.input                  = "infile:file"
pBedGetfasta.output                 = "outfile:file:{{in.infile | fn}}.fa"
pBedGetfasta.args.samtools          = params.samtools.value
pBedGetfasta.args.bedtools          = params.bedtools.value
pBedGetfasta.args.params            = Box({'name': True})
pBedGetfasta.args.ref               = params.ref.value
pBedGetfasta.tplenvs.runcmd         = runcmd.py
pBedGetfasta.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pBedGetfasta.beforeCmd              = checkref.fa.bash + buildref.fai.bash
pBedGetfasta.lang                   = params.python.value
pBedGetfasta.script                 = "file:scripts/bedtools/pBedGetfasta.py"


"""
@name:
	pBedClosest
@description:
	Similar to intersect, closest searches for overlapping features in A and B. In the event that no feature in B overlaps the current feature in A, closest will report the nearest (that is, least genomic distance from the start or end of A) feature in B. For example, one might want to find which is the closest gene to a significant GWAS polymorphism. Note that closest will report an overlapping feature as the closest that is, it does not restrict to closest non-overlapping feature. The following iconic cheatsheet summarizes the funcitonality available through the various optyions provided by the closest tool.
@input:
	`afile:file`:   The -a file
	`bfiles:files`: The -b files
@output:
	`outfile:file`: The result file
@args:
	`bin`:     The bedtools executable, default: "bedtools"
	`params`:  Other parameters for `bedtools closest`, default: ""
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedClosest = Proc()
pBedClosest.input  = "afile:file, bfiles:files"
pBedClosest.output = "outfile:file:{{afile | fn}}.bt"
pBedClosest.args   = { "bin": "bedtools", "params": "" }
pBedClosest.script = """
{{args.bin}} closest {{args.params}} -a "{{afile}}" -b {{bfiles | asquote}} > "{{outfile}}"
"""

"""
@name:
	pBedFlank
@description:
	`bedtools flank` will create two new flanking intervals for each interval in a BED/GFF/VCF file. Note that flank will restrict the created flanking intervals to the size of the chromosome (i.e. no start < 0 and no end > chromosome size).
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
pBedFlank = Proc()
pBedFlank.input  = "infile:file, gfile:file"
pBedFlank.output = "outfile:file:{{infile | fn}}.flank.bed"
pBedFlank.args   = { "bin": "bedtools", "params": "" }
pBedFlank.script = """
{{args.bin}} flank {{args.params}} -i "{{infile}}" -g "{{gfile}}" > "{{outfile}}"
"""

"""
@name:
	pBedIntersect
@description:
	By far, the most common question asked of two sets of genomic features is whether or not any of the features in the two sets overlap with one another. This is known as feature intersection. bedtools intersect allows one to screen for overlaps between two sets of genomic features. Moreover, it allows one to have fine control as to how the intersections are reported. bedtools intersect works with both BED/GFF/VCF and BAM files as input.
@input:
	`afile:file`:   The a file
	`bfiles:files`: The b files
@output:
	`outfile:file`: The result file
@args:
	`bin`:     The bedtools executable, default: "bedtools"
	`params`:  Other parameters for `bedtools intersect`, default: ""
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedIntersect = Proc()
pBedIntersect.input  = "afile:file, bfiles:files"
pBedIntersect.output = "outfile:file:{{afile|fn}}.intersect.bt"
pBedIntersect.args   = { "bin": "bedtools", "params": "" }
pBedIntersect.script = """
{{args.bin}} intersect -nonamecheck {{args.params}} -a "{{afile}}" -b {{bfiles | asquote}} > "{{outfile}}"
"""

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
pBedMakewindows = Proc()
pBedMakewindows.input  = "infile:file"
pBedMakewindows.output = "outfile:file:{{infile | fn}}.window.bed"
pBedMakewindows.args   = { "bin": "bedtools", "params": "" }
pBedMakewindows.script = """
{{args.bin}} makewindows {{args.params}} {{args.informat | lambda x: "-b" if x=="bed" else "-g"}} "{{infile}}" > "{{outfile}}"
"""

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
	`bin`:     The bedtools executable, default: "bedtools"
	`params`:  Other parameters for `bedtools merge`, default: ""
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedMerge = Proc()
pBedMerge.input  = "infile:file"
pBedMerge.output = "outfile:file:{{infile | fn}}.merged.bed"
pBedMerge.args   = { "bin": "bedtools", "params": "" }
pBedMerge.script = """
{{args.bin}} merge {{args.params}} -i "{{infile}}" > "{{outfile}}"
"""

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
pBedMultiinter = Proc()
pBedMultiinter.input  = "infiles:files"
pBedMultiinter.output = "outfile:file:{{infiles | [0] | fn}}.multiinter.bt"
pBedMultiinter.args   = { "bin": "bedtools", "params": "" }
pBedMultiinter.script = """
{{args.bin}} multiinter {{args.params}} -i {{infiles | asquote}} > "{{outfile}}"
"""

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
	`bin`:     The bedtools executable, default: "bedtools"
	`params`:  Other parameters for `bedtools random`, default: ""
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedRandom = Proc()
pBedRandom.input  = "gfile:file"
pBedRandom.output = "outfile:file:{{gfile | fn}}.random.bed"
pBedRandom.args   = { "bin": "bedtools", "params": "" }
pBedRandom.script = """
{{args.bin}} random {{args.params}} -g "{{gfile}}" > "{{outfile}}"
"""

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
pBedShift = Proc()
pBedShift.input  = "infile:file, gfile:file"
pBedShift.output = "outfile:file:{{infile | fn}}.shifted.bed"
pBedShift.args   = { "bin": "bedtools", "params": "" }
pBedShift.script = """
{{args.bin}} shift {{args.params}} -i "{{infile}}" -g "{{gfile}}" > "{{outfile}}"
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
pBedShuffle = Proc()
pBedShuffle.input  = "infile:file, gfile:file"
pBedShuffle.output = "outfile:file:{{infile | fn}}.shuffled.bed"
pBedShuffle.args   = { "bin": "bedtools", "params": "" }
pBedShuffle.script = """
{{args.bin}} shuffle {{args.params}} -i "{{infile}}" -g "{{gfile}}" > "{{outfile}}"
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
pBedSubtract = Proc()
pBedSubtract.input  = "afile:file, bfile:file"
pBedSubtract.output = "outfile:file:{{afile | fn}}.subtracted.bed"
pBedSubtract.args   = { "bin": "bedtools", "params": "" }
pBedSubtract.script = """
{{args.bin}} subtract {{args.params}} -a "{{afile}}" -b "{{bfile}}" > {{outfile}}
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
pBedWindow = Proc()
pBedWindow.input  = "afile:file, bfile:file"
pBedWindow.output = "outfile:file:{{afile | fn}}.window.bed"
pBedWindow.args   = { "bin": "bedtools", "params": "" }
pBedWindow.script = """
{{args.bin}} window {{args.params}} -a "{{afile}}" -b "{{bfile}}" > "{{outfile}}"
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
	`bin`:     The bedtools executable, default: "bedtools"
	`params`:  Other parameters for `bedtools genomecov`, default: "-bg"
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedGenomecov = Proc()
pBedGenomecov.input  = "infile:file"
pBedGenomecov.output = "outfile:file:{{infile | fn}}.genomecov.bt"
pBedGenomecov.args   = { "bin": "bedtools", "params": "-bg" }
pBedGenomecov.script = """
{{args.bin}} genomecov {{args.params}} -ibam "{{infile}}" > "{{outfile}}"
"""