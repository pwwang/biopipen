from pyppl import proc

#############################
# bedtools utilities        #
#############################

"""
@name:
	pGetfasta
@description:
	`bedtools getfasta` extracts sequences from a FASTA file for each of the intervals defined in a BED/GFF/VCF file.
@input:
	`infile:file`: The input bed file
	`fafile:file`: The input fasta file
@brings:
	`fafile`: "{{fafile | fn}}.fa*i", The fasta index file
@output:
	`outfile:file`: The generated fasta file
@args:
	`bin`:     The bedtools executable, default: "bedtools"
	`params`:  Other parameters for `bedtools getfasta`, default: ""
@requires:
	[bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pGetfasta = proc()
pGetfasta.input  = "infile:file, fafile:file"
pGetfasta.brings = {"fafile": "{{fafile | fn}}.fa*i"}
pGetfasta.output = "outfile:file:{{infile | fn}}.fa"
pGetfasta.args   = { "bin": "bedtools", "params": "-name" }
pGetfasta.script = """
{{proc.args.bin}} getfasta {{proc.args.params}} -fi "{{fafile}}" -bed "{{infile}}" > "{{outfile}}"
"""

"""
@name:
	pClosest
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
pClosest = proc()
pClosest.input  = "afile:file, bfiles:files"
pClosest.output = "outfile:file:{{afile | fn}}.bt"
pClosest.args   = { "bin": "bedtools", "params": "" }
pClosest.script = """
{{proc.args.bin}} closest {{proc.args.params}} -a "{{afile}}" -b {{bfiles | asquote}} > "{{outfile}}"
"""

"""
@name:
	pFlank
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
pFlank = proc()
pFlank.input  = "infile:file, gfile:file"
pFlank.output = "outfile:file:{{infile | fn}}.flank.bed"
pFlank.args   = { "bin": "bedtools", "params": "" }
pFlank.script = """
{{proc.args.bin}} flank {{proc.args.params}} -i "{{infile}}" -g "{{gfile}}" > "{{outfile}}"
"""

"""
@name:
	pIntersect
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
pIntersect = proc()
pIntersect.input  = "afile:file, bfiles:files"
pIntersect.output = "outfile:file:{{afile|fn}}.intersect.bt"
pIntersect.args   = { "bin": "bedtools", "params": "" }
pIntersect.script = """
{{proc.args.bin}} intersect -nonamecheck {{proc.args.params}} -a "{{afile}}" -b {{bfiles | asquote}} > "{{outfile}}"
"""

"""
@name:
	pMakewindows
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
pMakewindows = proc()
pMakewindows.input  = "infile:file"
pMakewindows.output = "outfile:file:{{infile | fn}}.window.bed"
pMakewindows.args   = { "bin": "bedtools", "params": "" }
pMakewindows.script = """
{{proc.args.bin}} makewindows {{proc.args.params}} {{proc.args.informat | lambda x: "-b" if x=="bed" else "-g"}} "{{infile}}" > "{{outfile}}"
"""

"""
@name:
	pMerge
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
pMerge = proc()
pMerge.input  = "infile:file"
pMerge.output = "outfile:file:{{infile | fn}}.merged.bed"
pMerge.args   = { "bin": "bedtools", "params": "" }
pMerge.script = """
{{proc.args.bin}} merge {{proc.args.params}} -i "{{infile}}" > "{{outfile}}"
"""

"""
@name:
	pMultiinter
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
pMultiinter = proc()
pMultiinter.input  = "infiles:files"
pMultiinter.output = "outfile:file:{{infiles | [0] | fn}}.multiinter.bt"
pMultiinter.args   = { "bin": "bedtools", "params": "" }
pMultiinter.script = """
{{proc.args.bin}} multiinter {{proc.args.params}} -i {{infiles | asquote}} > "{{outfile}}"
"""

"""
@name:
	pRandom
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
pRandom = proc()
pRandom.input  = "gfile:file"
pRandom.output = "outfile:file:{{gfile | fn}}.random.bed"
pRandom.args   = { "bin": "bedtools", "params": "" }
pRandom.script = """
{{proc.args.bin}} random {{proc.args.params}} -g "{{gfile}}" > "{{outfile}}"
"""

"""
@name:
	pShift
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
pShift = proc()
pShift.input  = "infile:file, gfile:file"
pShift.output = "outfile:file:{{infile | fn}}.shifted.bed"
pShift.args   = { "bin": "bedtools", "params": "" }
pShift.script = """
{{proc.args.bin}} shift {{proc.args.params}} -i "{{infile}}" -g "{{gfile}}" > "{{outfile}}"
"""

"""
@name:
	pShuffle
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
pShuffle = proc()
pShuffle.input  = "infile:file, gfile:file"
pShuffle.output = "outfile:file:{{infile | fn}}.shuffled.bed"
pShuffle.args   = { "bin": "bedtools", "params": "" }
pShuffle.script = """
{{proc.args.bin}} shuffle {{proc.args.params}} -i "{{infile}}" -g "{{gfile}}" > "{{outfile}}"
"""

"""
@name:
	pSubtract
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
pSubtract = proc()
pSubtract.input  = "afile:file, bfile:file"
pSubtract.output = "outfile:file:{{afile | fn}}.subtracted.bed"
pSubtract.args   = { "bin": "bedtools", "params": "" }
pSubtract.script = """
{{proc.args.bin}} subtract {{proc.args.params}} -a "{{afile}}" -b "{{bfile}}" > {{outfile}}
"""

"""
@name:
	pWindow
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
pWindow = proc()
pWindow.input  = "afile:file, bfile:file"
pWindow.output = "outfile:file:{{afile | fn}}.window.bed"
pWindow.args   = { "bin": "bedtools", "params": "" }
pWindow.script = """
{{proc.args.bin}} window {{proc.args.params}} -a "{{afile}}" -b "{{bfile}}" > "{{outfile}}"
"""

"""
@name:
	pGenomecov
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
pGenomecov = proc()
pGenomecov.input  = "infile:file"
pGenomecov.output = "outfile:file:{{infile | fn}}.genomecov.bt"
pGenomecov.args   = { "bin": "bedtools", "params": "-bg" }
pGenomecov.script = """
{{proc.args.bin}} genomecov {{proc.args.params}} -ibam "{{infile}}" > "{{outfile}}"
"""