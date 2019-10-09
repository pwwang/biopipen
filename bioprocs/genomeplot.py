"""Plot genomic features using Gviz R package"""
# https://bioconductor.org/packages/devel/bioc/vignettes/Gviz/inst/doc/Gviz.pdf
from pyppl import Proc, Box
from . import params, rimport
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())


@procfactory
def _pInteractionTrack():
	"""
	@name:
		pInteractionTrack
	@description:
		Gererate genomic interaction track for Gviz
	@input:
		`name`: The name of the track
		`infile:file`: The input file.
			- See the `type` argument for `makeGenomicInteractionsFromFile` from `GenomicInteractions` r-package
		`region`: the region, just chromosome!
	@output:
		`outfile:file`: The dumped track data
	@args:
		`intype`: Input file type. Default: auto
			- Identified by extension
			- One of "chiapet.tool", "bed12", "bedpe", "hiclib", "homer", "bam", "two.bams".
		`params`: The display params
	"""
	pInteractionTrack             = Proc(desc = 'Gererate genomic interaction track for Gviz.')
	pInteractionTrack.input       = "name, infile:file, region"
	pInteractionTrack.output      = "outfile:file:interactionTrack_{{i.name}}_{{i.region | lambda x: x.replace(':', '-')}}.gviztrack"
	pInteractionTrack.args.intype = "auto"
	pInteractionTrack.args.params = Box({
		'background.title': "#333333",
		'col'             : 'NULL',
		'col.outside'     : "#99beff",
		'anchor.height'   : 0.06,
		'col.anchors.fill': '#99beff',
		'col.anchors.line': '#AAAAAA',
		'col.interactions': '#f98d85',
		'size'            : 2
	})
	pInteractionTrack.lang        = params.Rscript.value
	pInteractionTrack.script      = "file:scripts/genomeplot/pInteractionTrack.r"
	return pInteractionTrack

@procfactory
def _pGeneTrack():
	"""
	@name:
		pGeneTrack
	@description:
		Generate gene track using ucsc data source
	@input:
		`name`:   The name of the track
	@output:
		`outfile:file`: The file to save the track
	@args:
		`genome`: The genome
		`params`: use `displayPars(UcscTrack(genome="mm9", chromosome="chrM", track="knownGene"))` to see all available args
	@requires:
		[r-Gviz](https://rdrr.io/bioc/Gviz)
	"""
	pGeneTrack             = Proc(desc = 'Generate gene track data for pGenomePlot.')
	pGeneTrack.input       = "name, region"
	pGeneTrack.output      = "outfile:file:geneTrack_{{i.name}}_{{i.region | lambda x: x.replace(':', '-')}}.gviztrack"
	pGeneTrack.args.genome = params.genome.value
	pGeneTrack.args.params = Box()
	pGeneTrack.lang        = params.Rscript.value
	pGeneTrack.script      = "file:scripts/genomeplot/pGeneTrack.r"
	return pGeneTrack

@procfactory
def _pAnnoTrack():
	"""
	@name:
		pAnnoTrack
	@description:
		The annotation track of Gviz
	@input:
		`name`:        the name of the track
		`infile:file`: the file for the track. (wig, bigWig or bedGraph, bam, need to be indexed!)
		`chrom`:       the chrom
	@output:
		`outfile:file`:the rds file for the track
	@args:
		`genome`: The genome
		`params`:  See `displayPars(DataTrack())` for all available display params
	@requires:
		[r-Gviz](https://rdrr.io/bioc/Gviz/man/DataTrack-class.html)
	"""
	pAnnoTrack             = Proc(desc = 'Generate annotation track for pGenomePlot.')
	pAnnoTrack.input       = "name, infile:file, chrom"
	pAnnoTrack.output      = "outfile:file:dataTrack_{{i.name}}_{{i.chrom}}.gviztrack"
	pAnnoTrack.args.genome = params.genome.value
	pAnnoTrack.args.params = Box()
	pAnnoTrack.lang        = params.Rscript.value
	pAnnoTrack.script      = "file:scripts/genomeplot/pAnnoTrack.r"
	return pAnnoTrack

@procfactory
def _pDataTrack():
	"""
	@name:
		pDataTrack
	@description:
		The data track of Gviz
	@input:
		`name`:        the name of the track
		`infile:file`: the file for the track. (wig, bigWig or bedGraph, bam, need to be indexed!)
		`chrom`:       the chrom
	@output:
		`outfile:file`:the rds file for the track
	@args:
		`genome`: The genome
		`params`:  See `displayPars(DataTrack())` for all available display params
	@requires:
		[r-Gviz](https://rdrr.io/bioc/Gviz/man/DataTrack-class.html)
	"""
	pDataTrack             = Proc(desc = 'Generate data track for pGenomePlot.')
	pDataTrack.input       = "name, infile:file, chrom"
	pDataTrack.output      = "outfile:file:dataTrack_{{i.name}}_{{i.chrom}}.gviztrack"
	pDataTrack.args.genome = params.genome.value
	pDataTrack.args.params = Box()
	pDataTrack.lang        = params.Rscript.value
	pDataTrack.script      = "file:scripts/genomeplot/pDataTrack.r"
	return pDataTrack

@procfactory
def _pUcscTrack():
	"""
	@name:
		pUcscTrack
	@description:
		Generate track from ucsc
	@input:
		`name`     : the name of the track
		`track`    : the UCSC track
		`trackType`: the Gviz track
		`region`   : the region
	@output:
		`outfile:file`:the dumped track
	@args:
		`params`: use `displayPars(UcscTrack(genome="mm9", chromosome="chrM", track="knownGene"))` to see all available args.
	@requires:
		[r-Gviz](https://rdrr.io/bioc/Gviz)
	"""
	pUcscTrack             = Proc(desc = "Generate track from UCSC data.")
	pUcscTrack.input       = "name, track, trackType, region"
	pUcscTrack.output      = "outfile:file:ucscTrack_{{i.name}}_{{i.region | lambda x: x.replace(':', '-')}}.gviztrack"
	pUcscTrack.lang        = params.Rscript.value
	pUcscTrack.args.genome = params.genome.value
	pUcscTrack.args.params = Box()
	pUcscTrack.script      = "file:scripts/genomeplot/pUcscTrack.r"
	return pUcscTrack

@procfactory
def _pGenomePlot():
	"""
	@name:
		pGenomePlot
	@description:
		plot the genomic features
	@input:
		`trkfiles:files`: the list of track dumped files
		`region`:         the region, in format of `chr1:1-1000`
		`highlight`:      the highlight regions, informat of start1-end1; start2-end2; ...
	@output:
		`outfile:file`:   the figure
	@args:
		`genome`  : The genome
		`showIdeo`: Show ideogram track? Default: True
		`showAxis`: Show axis? Default: True
		`showGenes`: Show geneTrack? Default: True
		`params`:   The params
			- `genneral`:  General params for plotTracks
			- `geneTrack`: The params for geneTrack
	@requires:
		[r-Gviz](https://rdrr.io/bioc/Gviz)
	"""
	pGenomePlot                = Proc(desc = 'Plot genome elements.')
	pGenomePlot.input          = "trkfiles:files, region, highlight"
	pGenomePlot.output         = "outfile:file:genomeplot_{{i.region | lambda x: x.replace(':', '-')}}.png"
	pGenomePlot.args.genome    = params.genome.value
	pGenomePlot.args.ideoTrack = params.cytoband.value # a file from ucsc (http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz); or True to download it in runtime; or False to disable it
	pGenomePlot.args.axisTrack = True
	pGenomePlot.args.geneTrack = params.refgene.value # or False to disable
	pGenomePlot.args.params    = Box({
		'title.width': .8
	})
	pGenomePlot.args.scheme  = Box(
		GenomeAxisTrack = Box({
			'fontsize': 8
		}),
		GeneRegionTrack = Box({
			'fontsize'         : 12,
			'fill'             : "salmon",
			'col'              : "NULL",
			'showId'           : True,
			'shape'            : ["smallArrow", "arrow"],
			'size'             : .5,
			'background.title' : "#333333",
			"arrowHeadWidth"   : 30,
			"arrowHeadMaxWidth": 40
		}),
		AnnotationTrack = Box({
			'fontsize'        : 12,
			'background.title': "#333333",
			'fontcolor.item'  : "#333333"
		})
	)
	pGenomePlot.args.devpars = Box(res = 300, height = 300, width = 2000)
	pGenomePlot.envs.rimport = rimport
	pGenomePlot.lang         = params.Rscript.value
	pGenomePlot.script       = "file:scripts/genomeplot/pGenomePlot.r"
	return pGenomePlot


