from pyppl import Proc, Box
from . import params
"""
Plot genomic features using Gviz R package
https://bioconductor.org/packages/devel/bioc/vignettes/Gviz/inst/doc/Gviz.pdf
"""

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
pInteractionTrack.output      = "outfile:file:interactionTrack_{{in.name}}_{{in.region | lambda x: x.replace(':', '-')}}.gviztrack"
pInteractionTrack.args.intype = "auto"
pInteractionTrack.args.params = Box({
    'background.title': "brown",
    'col'             : 'NULL',
    'fill'            : '#FFD58A',
    'col.outside'     : "lightblue",
    'anchor.height'   : 0.06,
    'col.anchors.fill': 'lightblue',
	'col.anchors.line': '#AAAAAA'
})
pInteractionTrack.lang        = params.Rscript.value
pInteractionTrack.script      = "file:scripts/genomeplot/pInteractionTrack.r"

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
pGeneTrack.output      = "outfile:file:geneTrack_{{in.name}}_{{in.region | lambda x: x.replace(':', '-')}}.gviztrack"
pGeneTrack.args.genome = params.genome.value
pGeneTrack.args.params = Box({
    "arrowHeadWidth"   : 30,
    "arrowHeadMaxWidth": 40,
    "shape"            : "arrow",
    "showId"           : True,
    'background.title' : "brown",
    'col'              : 'NULL',
    'fill'             : '#FFD58A'
})
pGeneTrack.lang   = params.Rscript.value
pGeneTrack.script = "file:scripts/genomeplot/pGeneTrack.r"

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
pAnnoTrack.brings      = {"infile": ["{{in.infile | fn}}.bai", "{{in.infile | bn}}.bai", "{{in.infile | bn}}"]}
pAnnoTrack.output      = "outfile:file:dataTrack_{{in.name}}_{{in.chrom}}.gviztrack"
pAnnoTrack.args.genome = params.genome.value
pAnnoTrack.args.params = Box({
    'background.title': "brown",
    'col'             : 'NULL',
    'fill'            : '#FFD58A'
})
pAnnoTrack.lang        = params.Rscript.value
pAnnoTrack.script      = "file:scripts/genomeplot/pAnnoTrack.r"

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
pDataTrack.brings      = {"infile": ["{{in.infile | fn}}.bai", "{{in.infile | bn}}.bai", "{{in.infile | bn}}"]}
pDataTrack.output      = "outfile:file:dataTrack_{{in.name}}_{{in.chrom}}.gviztrack"
pDataTrack.args.genome = params.genome.value
pDataTrack.args.params = Box({
    'background.title': "brown",
    'col'             : 'NULL',
    'fill'            : '#FFD58A'
})
pDataTrack.lang        = params.Rscript.value
pDataTrack.script      = "file:scripts/genomeplot/pDataTrack.r"

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
    use `displayPars(UcscTrack(genome="mm9", chromosome="chrM", track="knownGene"))` to see all available args.
@requires:
    [r-Gviz](https://rdrr.io/bioc/Gviz)
"""
pUcscTrack             = Proc(desc = "Generate track from UCSC data.")
pUcscTrack.input       = "name, track, trackType, region"
pUcscTrack.output      = "outfile:file:ucscTrack_{{in.name}}_{{in.region | lambda x: x.replace(':', '-')}}.gviztrack"
pUcscTrack.lang        = params.Rscript.value
pUcscTrack.args.genome = params.genome.value
pUcscTrack.args.params = Box({
    'background.title': "brown",
    'col'             : 'NULL',
    'fill'            : '#FFD58A'
})
pUcscTrack.script      = "file:scripts/genomeplot/pUcscTrack.r"

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
pGenomePlot.output         = "outfile:file:genomeplot_{{in.region | lambda x: x.replace(':', '-')}}.png"
pGenomePlot.args.genome    = params.genome.value
pGenomePlot.args.ideoTrack = params.cytoband.value # a file from ucsc (http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz) or True to download it in runtime
pGenomePlot.args.axisTrack = True
pGenomePlot.args.geneTrack = params.refgene.value # or false to disable
pGenomePlot.args.params    = Box(
    general   = Box(),
    axisTrack = Box(

    ),
    ideoTrack = Box(

    ),
    geneTrack = Box({
        "shape"            : "arrow",
        "showId"           : True,
        'background.title' : "brown",
        'col'              : 'NULL',
        'fill'             : '#FFD58A'
    })
)
pGenomePlot.args.devpars   = Box(res = 300, height = 1000, width = 4000)
pGenomePlot.lang           = params.Rscript.value
pGenomePlot.script         = "file:scripts/genomeplot/pGenomePlot.r"


