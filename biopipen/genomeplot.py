"""Plot genomic features using Gviz R package"""
# https://bioconductor.org/packages/devel/bioc/vignettes/Gviz/inst/doc/Gviz.pdf
from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name

pInteractionTrack = proc_factory(
    desc='Gererate genomic interaction track for Gviz.',
    config=Diot(annotate="""
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
    """),
    input="name, infile:file, region",
    output="""outfile:file:interactionTrack_{{i.name}}_{{
        i.region | .replace: ':', '-'}}.gviztrack""",
    lang=opts.Rscript,
    args=Diot(
        intype="auto",
        params=Diot({
            'background.title': "#333333",
            'col': 'NULL',
            'col.outside': "#99beff",
            'anchor.height': 0.06,
            'col.anchors.fill': '#99beff',
            'col.anchors.line': '#AAAAAA',
            'col.interactions': '#f98d85',
            'size': 2
        })
    )
)

pGeneTrack = proc_factory(
    desc='Generate gene track data for pGenomePlot.',
    config=Diot(annotate="""
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
    """),
    input="name, region",
    output="""outfile:file:geneTrack_{{i.name}}_{{
        i.region | .replace: ':', '-'}}.gviztrack""",
    lang=opts.Rscript,
    args=Diot(
        genome=opts.genome,
        params=Diot(),
    )
)

pAnnoTrack = proc_factory(
    desc='Generate annotation track for pGenomePlot.',
    config=Diot(annotate="""
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
    """),
    input="name, infile:file, chrom",
    output="outfile:file:dataTrack_{{i.name}}_{{i.chrom}}.gviztrack",
    lang=opts.Rscript,
    args=Diot(
        genome=opts.genome,
        params=Diot(),
    )
)

pDataTrack = proc_factory(
    desc='Generate data track for pGenomePlot.',
    config=Diot(annotate="""
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
    """),
    input="name, infile:file, chrom",
    output="outfile:file:dataTrack_{{i.name}}_{{i.chrom}}.gviztrack",
    lang=opts.Rscript,
    args=Diot(
        genome=opts.genome,
        params=Diot(),
    )
)

pUcscTrack = proc_factory(
    desc="Generate track from UCSC data.",
    config=Diot(annotate="""
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
    """),
    input="name, track, trackType, region",
    output="""outfile:file:ucscTrack_{{i.name}}_{{
        i.region | .replace: ':', '-'}}.gviztrack""",
    lang=opts.Rscript,
    args=Diot(
        genome=opts.genome,
        params=Diot(),
    )
)

pGenomePlot = proc_factory(
    desc='Plot genome elements.',
    config=Diot(annotate="""
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
    """),
    input="trkfiles:files, region, highlight",
    output="outfile:file:genomeplot_{{i.region | .replace: ':', '-'}}.png",
    lang=opts.Rscript,
    args=Diot(
        genome=opts.genome,
        ideoTrack=opts.cytoband,
        # a file from ucsc
        # http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
        # or True to download it in runtime; or False to disable it,
        axisTrack=True,
        geneTrack=opts.refgene,  # or False to disable,
        params=Diot({'title.width': .8}),
        scheme=Diot(GenomeAxisTrack=Diot({'fontsize': 8}),
                    GeneRegionTrack=Diot({
                        'fontsize': 12,
                        'fill': "salmon",
                        'col': "NULL",
                        'showId': True,
                        'shape': ["smallArrow", "arrow"],
                        'size': .5,
                        'background.title': "#333333",
                        "arrowHeadWidth": 30,
                        "arrowHeadMaxWidth": 40
                    }),
                    AnnotationTrack=Diot({
                        'fontsize': 12,
                        'background.title': "#333333",
                        'fontcolor.item': "#333333"
                    })),
        devpars=Diot(res=300, height=300, width=2000),
    ),
)
