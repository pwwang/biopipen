from pyppl import proc

"""
Plot genomic features using Gviz R package
https://bioconductor.org/packages/devel/bioc/vignettes/Gviz/inst/doc/Gviz.pdf
"""

"""
@name:
    pGeneTrack
@description:
    Generate the gene track using ucsc data source
@input:
    `name`:   The name of the track
    `genome`: The genome
    `chrom`:  The chromosome
    `from`:   The start
    `to`:     The end
@output:
    `outfile:file`: The file to save the track
@args:
    use `displayPars(UcscTrack(genome="mm9", chromosome="chrM", track="knownGene"))` to see all available args
@requires:
    [r-Gviz](https://rdrr.io/bioc/Gviz)
"""
pGeneTrack = proc ()
pGeneTrack.input   = "name, genome, chrom, from, to"
pGeneTrack.output  = "outfile:file:geneTrack.{{name}}.{{#}}.rds, gout:{{genome}}, cout:{{chrom}}, fout:{{from}}, tout: {{to}}"
pGeneTrack.lang    = "Rscript"
pGeneTrack.args    = {
    "arrowHeadWidth":          "30",
    "arrowHeadMaxWidth":       "40",
    "cex.group":               "0.6",
    "cex":                     "1",
    "col.line":                "darkgray",
    "col":                     "darkgray",
    "featureAnnotation":       "NULL",
    "fill":                    "lightblue",
    "fontfamily.group":        "sans",
    "fontcolor.group":         "#808080",
    "fontcolor.item":          "white",
    "fontface.group":          "2",
    "fontsize.group":          "12",
    "groupAnnotation":         "NULL",
    "just.group":              "left",
    "lex":                     "1",
    "lineheight":              "1",
    "lty":                     "solid",
    "lwd.baseline":            "NULL",
    "lwd":                     "1",
    "mergeGroups":             "FALSE",
    "min.height":              "3",
    "min.width":               "1",
    "rotation":                "0",
    "rotation.group":          "0",
    "rotation.item":           "0",
    "shape":                   "arrow",
    "showFeatureId":           "NULL",
    "showId":                  "TRUE",
    "showOverplotting":        "FALSE",
    "size":                    "1",
    "stackHeight":             "0.75",
    "reverseStacking":         "FALSE",
    "alpha":                   "1",
    "alpha.title":             "NULL",
    "background.panel":        "transparent",
    "background.title":        "lightgray",
    "cex.axis":                "NULL",
    "cex.title":               "NULL",
    "col.axis":                "white",
    "col.border.title":        "white",
    "col.frame":               "lightgray",
    "col.grid":                "#808080",
    "col.symbol":              "NULL",
    "fontcolor.title":         "white",
    "collapse":                "TRUE",
    "fontcolor":               "black",
    "fontface.title":          "2",
    "fontface":                "1",
    "fontfamily.title":        "sans",
    "fontfamily":              "sans",
    "fontsize":                "12",
    "frame":                   "FALSE",
    "grid":                    "FALSE",
    "h":                       "-1",
    "lty.grid":                "solid",
    "lwd.border.title":        "1",
    "lwd.title":               "1",
    "lwd.grid":                "1",
    "min.distance":            "1",
    "reverseStrand":           "FALSE",
    "rotation.title":          "90",
    "showAxis":                "TRUE",
    "showTitle":               "TRUE",
    "v":                       "-1",
    "lwd.border":              "1"
}
pGeneTrack.script  = """
library (Gviz)
geneTrack = UcscTrack (
    genome          = "{{genome}}",
    chromosome      = "{{chrom}}",
    from            = {{from}},
    to              = {{to}},
    name            = "{{name}}",
    trackType       = "GeneRegionTrack",
    table           = "refGene",
    track           = "refGene",
    rstarts         = "exonStarts",
    rends           = "exonEnds",
    gene            = "name",
    symbol          = "name2",
    transcript      = "name",
    strand          = "strand",
{{ proc.args | \
  (lambda x: ",\\n".join( \
     [ \
        "\\t%s\\t= %s" % (key, ('"'+ val +'"' if key in [ \
            "col.line", "col", "fill", "fontfamily.group", "fontcolor.group", "fontcolor.item", \
            "just.group", "lty", "shape", "background.panel", "background.title", "col.border.title", \
            "col.frame", "col.grid", "fontcolor.title", "fontcolor", "fontfamily.title", "fontfamily", \
            "lty.grid", "col.axis"] and val != "NULL" else val))  \
        for key, val in x.iteritems() \
     ] \
  ))(_) }}
)
saveRDS (geneTrack, "{{outfile}}")
"""

"""
@name:
    pAnnoTrack
@description:
    Generate annotation track
@input:
    `infile:file`: the file for the track
    `name`:        the name of the track
    `genome`:      the genome
    `chrom`:       the chromosome
    `from`:        the start position to display
    `to`:          the end position to display
@output:
    `outfile:file`:the dumped track
@args:
    use `displayPars(AnnotationTrack())` to see all available args.
@requires:
    [r-Gviz](https://rdrr.io/bioc/Gviz)
"""
pAnnoTrack = proc ()
pAnnoTrack.input   = "infile:file, name, genome, chrom, from, to"
pAnnoTrack.output  = "outfile:file:annoTrack.{{name}}.{{#}}.rds, gout:{{genome}}, cout:{{chrom}}, fout:{{from}}, tout: {{to}}"
pAnnoTrack.lang  = "Rscript"
pAnnoTrack.args    = {
    "group":                   "NULL",
    "arrowHeadWidth":          "30",
    "arrowHeadMaxWidth":       "40",
    "cex.group":               "0.6",
    "cex":                     "1",
    "col.line":                "darkgray",
    "col":                     "darkgray",
    "featureAnnotation":       "id",
    "fill":                    "lightblue",
    "fontfamily.group":        "sans",
    "fontcolor.group":         "#808080",
    "fontcolor.item":          "#555555",
    "fontface.group":          "2",
    "fontsize.group":          "12",
    "groupAnnotation":         "NULL",
    "just.group":              "left",
    "lex":                     "1",
    "lineheight":              "1",
    "lty":                     "solid",
    "lwd.baseline":            "NULL",
    "lwd":                     "1",
    "mergeGroups":             "FALSE",
    "min.height":              "3",
    "min.width":               "1",
    "rotation":                "0",
    "rotation.group":          "0",
    "rotation.item":           "0",
    "shape":                   "arrow",
    "showFeatureId":           "NULL",
    "showId":                  "FALSE",
    "showOverplotting":        "FALSE",
    "size":                    "1",
    "stackHeight":             "0.75",
    "reverseStacking":         "FALSE",
    "alpha":                   "1",
    "alpha.title":             "NULL",
    "background.panel":        "transparent",
    "background.title":        "lightgray",
    "cex.axis":                "NULL",
    "cex.title":               "NULL",
    "col.axis":                "white",
    "col.border.title":        "white",
    "col.frame":               "lightgray",
    "col.grid":                "#808080",
    "col.symbol":              "NULL",
    "fontcolor.title":         "white",
    "collapse":                "TRUE",
    "fontcolor":               "black",
    "fontface.title":          "2",
    "fontface":                "1",
    "fontfamily.title":        "sans",
    "fontfamily":              "sans",
    "fontsize":                "12",
    "frame":                   "FALSE",
    "grid":                    "FALSE",
    "h":                       "-1",
    "lty.grid":                "solid",
    "lwd.border.title":        "1",
    "lwd.title":               "1",
    "lwd.grid":                "1",
    "min.distance":            "1",
    "reverseStrand":           "FALSE",
    "rotation.title":          "90",
    "showAxis":                "TRUE",
    "showTitle":               "TRUE",
    "v":                       "-1",
    "lwd.border":              "1"
}
pAnnoTrack.script  = """
library (Gviz)
annoTrack = AnnotationTrack(
    range = "{{infile}}",
    genome = "{{genome}}",
    name = "{{name}}",
    chromosome = "{{chrom}}",
{{ proc.args | \
  (lambda x: ",\\n".join( \
     [ \
        "\\t%s\\t= %s" % (key, ('"'+ val +'"' if key in [ \
            "col.line", "col", "fill", "fontfamily.group", "fontcolor.group", "fontcolor.item", \
            "just.group", "lty", "shape", "background.panel", "background.title", "col.border.title", \
            "col.frame", "col.grid", "fontcolor.title", "fontcolor", "fontfamily.title", "fontfamily", \
            "lty.grid", "col.axis", "featureAnnotation"] and val != "NULL" else val))  \
        for key, val in x.iteritems() \
     ] \
  ))(_) }}
)
# the group information
if (is.null({{proc.args.group}})) {
    group (annoTrack) = as.character (seq(length(annoTrack)))
} else {
    # a function to generate groups
    group (annoTrack) = {{proc.args.group}} (annoTrack)
}
saveRDS (annoTrack, "{{outfile}}")
"""

"""
@name:
    pDataTrack
@description:
    The data track of Gviz
@input:
    `infile:file`: the file for the track
    `name`:        the name of the track
    `genome`:      the genome
    `chrom`:       the chromosome
    `from`:        the start position to display
    `to`:          the end position to display
@output:
    `outfile:file`:the rds file for the track
    `gout`:        the genome
    `cout`:        the chromosome
    `fout`:        the start
    `tout`:        the end
@args:
    See `displayPars(DataTrack())` for all available display params
    Quote all params!
@requires:
    [r-Gviz](https://rdrr.io/bioc/Gviz/man/DataTrack-class.html)
"""
pDataTrack = proc ()
pDataTrack.input   = "infile:file, name, genome, chrom, from, to"
pDataTrack.output  = "outfile:file:dataTrack.{{name}}.{{#}}.rds, gout:{{genome}}, cout:{{chrom}}, fout:{{from}}, tout: {{to}}"
pDataTrack.args    = {
    "start":                   "NULL",
    "end":                     "NULL",
    "width":                   "NULL",
    "aggregateGroups":         "FALSE",
    "alpha.confint":           "0.3",
    "amount":                  "NULL",
    "baseline":                "NULL",
    "box.legend":              "FALSE",
    "box.ratio":               "1",
    "box.width":               "NULL",
    "grid":                    "FALSE",
    "cex.legend":              "0.8",
    "cex.sampleNames":         "NULL",
    "cex":                     "0.7",
    "coef":                    "1.5",
    "col.baseline":            "NULL",
    "col.confint":             "NA",
    "col.histogram":           "#808080",
    "col.horizon":             "NA",
    "col.mountain":            "NULL",
    "col.sampleNames":         "white",
    "col":                     "c(\"#0080ff\", \"#ff00ff\", \"darkgreen\", \"#ff0000\", \"orange\", \"#00ff00\", \"brown\")",
    "collapse":                "FALSE",
    "degree":                  "1",
    "do.out":                  "TRUE",
    "evaluation":              "50",
    "factor":                  "0.5",
    "family":                  "symmetric",
    "fill.confint":            "NULL",
    "fill.histogram":          "NULL",
    "fill.horizon":            "c(\"#B41414\", \"#E03231\", \"#F7A99C\", \"#9FC8DC\", \"#468CC8\", \"#0165B3\")",
    "fill.mountain":           "c(\"#CCFFFF\", \"#FFCCFF\")",
    "fontface.legend":         "NULL",
    "fontfamily.legend":       "NULL",
    "fontsize.legend":         "NULL",
    "fontcolor.legend":        "#808080",
    "gradient":                "c(\"#F7FBFF\", \"#DEEBF7\", \"#C6DBEF\", \"#9ECAE1\", \"#6BAED6\", \"#4292C6\", \"#2171B5\", \"#08519C\", \"#08306B\")",
    "groups":                  "NULL",
    "horizon.origin":          "0",
    "horizon.scale":           "NULL",
    "jitter.x":                "FALSE",
    "jitter.y":                "FALSE",
    "levels.fos":              "NULL",
    "legend":                  "TRUE",
    "lineheight.legend":       "NULL",
    "lty.baseline":            "NULL",
    "lty.mountain":            "NULL",
    "lwd.baseline":            "NULL",
    "lwd.mountain":            "NULL",
    "min.distance":            "0",
    "na.rm":                   "FALSE",
    "ncolor":                  "100",
    "notch.frac":              "0.5",
    "notch":                   "FALSE",
    "pch":                     "20",
    "separator":               "0",
    "showColorBar":            "TRUE",
    "showSampleNames":         "FALSE",
    "size":                    "NULL",
    "span":                    "0.2",
    "stackedBars":             "TRUE",
    "stats":                   """function (x, coef = 1.5, do.conf = TRUE, do.out = TRUE) \\n\
{                                                                                         \\n\
    if (coef < 0)                                                                         \\n\
        stop("'coef' must not be negative")                                               \\n\
    nna <- !is.na(x)                                                                      \\n\
    n <- sum(nna)                                                                         \\n\
    stats <- stats::fivenum(x, na.rm = TRUE)                                              \\n\
    iqr <- diff(stats[c(2, 4)])                                                           \\n\
    if (coef == 0)                                                                        \\n\
        do.out <- FALSE                                                                   \\n\
    else {                                                                                \\n\
        out <- if (!is.na(iqr)) {                                                         \\n\
            x < (stats[2L] - coef * iqr) | x > (stats[4L] + coef *                        \\n\
                iqr)                                                                      \\n\
        }                                                                                 \\n\
        else !is.finite(x)                                                                \\n\
        if (any(out[nna], na.rm = TRUE))                                                  \\n\
            stats[c(1, 5)] <- range(x[!out], na.rm = TRUE)                                \\n\
    }                                                                                     \\n\
    conf <- if (do.conf)                                                                  \\n\
        stats[3L] + c(-1.58, 1.58) * iqr/sqrt(n)                                          \\n\
    list(stats = stats, n = n, conf = conf, out = if (do.out) x[out &                     \\n\
        nna] else numeric())                                                              \\n\
}""",
    "transformation":          "NULL",
    "type":                    "p",
    "varwidth":                "FALSE",
    "window":                  "NULL",
    "windowSize":              "NULL",
    "ylim":                    "NULL",
    "alpha":                   "1",
    "alpha.title":             "NULL",
    "background.panel":        "transparent",
    "background.title":        "lightgray",
    "cex.axis":                "NULL",
    "cex.title":               "NULL",
    "col.axis":                "white",
    "col.border.title":        "white",
    "col.frame":               "lightgray",
    "col.grid":                "#808080",
    "col.line":                "NULL",
    "col.symbol":              "NULL",
    "fontcolor.title":         "white",
    "fill":                    "lightgray",
    "fontcolor":               "black",
    "fontface.title":          "2",
    "fontface":                "1",
    "fontfamily.title":        "sans",
    "fontfamily":              "sans",
    "fontsize":                "12",
    "frame":                   "FALSE",
    "h":                       "-1",
    "lineheight":              "1",
    "lty.grid":                "solid",
    "lty":                     "solid",
    "lwd.border.title":        "1",
    "lwd.title":               "1",
    "lwd.grid":                "1",
    "lwd":                     "1",
    "min.height":              "3",
    "min.width":               "1",
    "reverseStrand":           "FALSE",
    "rotation.title":          "90",
    "rotation":                "0",
    "showAxis":                "TRUE",
    "showTitle":               "TRUE",
    "v":                       "-1",
    "lwd.border":              "1"
}
pDataTrack.lang    = "Rscript"
pDataTrack.script  = """
library (Gviz)
dataTrack = DataTrack(
    range = "{{infile}}",
    genome = "{{genome}}",
    name = "{{name}}",
    chromosome = "{{chrom}}",
{{ proc.args | \
  (lambda x: ",\\n".join( \
     [ \
        "\\t%s\\t= %s" % (key, ('"'+ val +'"' if key in [ \
            "col.histogram", "col.sampleNames", "family", "fontcolor.legend", "type", \
            "background.panel", "background.title", "background.panel", "col.axis", \
            "col.border.title", "col.frame", "col.grid", "fontcolor.title", "fill", "fontcolor", \
            "fontfamily.title", "fontfamily", "lty.grid", "lty"] and val != "NULL" else val))  \
        for key, val in x.iteritems() \
     ] \
  ))(_) }}
)
saveRDS (dataTrack, "{{outfile}}")
"""

"""
@name:
    pUcscTrack
@description:
    Generate track from ucsc
@input:
    `ucscTrack`:   the track to fetch from ucsc. [Avialable tracks](http://genome.ucsc.edu/cgi-bin/hgTables?command=start)
    `table`:       the table from ucsc. [Available table](http://genome.ucsc.edu/cgi-bin/hgTables?command=start)
    `gvizTrack`:   the object track to generate. One of "AnnotationTrack", "GeneRegionTrack", "DataTrack", "GenomeAxisTrack"
    `name`:        the name of the track
    `genome`:      the genome
    `chrom`:       the chromosome
    `from`:        the start position to display
    `to`:          the end position to display
@output:
    `outfile:file`:the dumped track
@args:
    use `displayPars(UcscTrack(genome="mm9", chromosome="chrM", track="knownGene"))` to see all available args.
@requires:
    [r-Gviz](https://rdrr.io/bioc/Gviz)
"""
pUcscTrack = proc ()
pUcscTrack.input   = "ucscTrack, table, gvizTrack, name, genome, chrom, from, to"
pUcscTrack.output  = "outfile:file:ucscTrack.{{name}}.{{#}}.rds, gout:{{genome}}, cout:{{chrom}}, fout:{{from}}, tout: {{to}}"
pUcscTrack.lang    = "Rscript"
pUcscTrack.args    = {
    "arrowHeadWidth":          "30",
    "arrowHeadMaxWidth":       "40",
    "cex.group":               "0.6",
    "cex":                     "1",
    "col.line":                "darkgray",
    "col":                     "darkgray",
    "featureAnnotation":       "NULL",
    "fill":                    "lightblue",
    "fontfamily.group":        "sans",
    "fontcolor.group":         "#808080",
    "fontcolor.item":          "white",
    "fontface.group":          "2",
    "fontsize.group":          "12",
    "groupAnnotation":         "NULL",
    "just.group":              "left",
    "lex":                     "1",
    "lineheight":              "1",
    "lty":                     "solid",
    "lwd.baseline":            "NULL",
    "lwd":                     "1",
    "mergeGroups":             "FALSE",
    "min.height":              "3",
    "min.width":               "1",
    "rotation":                "0",
    "rotation.group":          "0",
    "rotation.item":           "0",
    "shape":                   "arrow",
    "showFeatureId":           "NULL",
    "showId":                  "NULL",
    "showOverplotting":        "FALSE",
    "size":                    "1",
    "stackHeight":             "0.75",
    "reverseStacking":         "FALSE",
    "alpha":                   "1",
    "alpha.title":             "NULL",
    "background.panel":        "transparent",
    "background.title":        "lightgray",
    "cex.axis":                "NULL",
    "cex.title":               "NULL",
    "col.axis":                "white",
    "col.border.title":        "white",
    "col.frame":               "lightgray",
    "col.grid":                "#808080",
    "col.symbol":              "NULL",
    "fontcolor.title":         "white",
    "collapse":                "TRUE",
    "fontcolor":               "black",
    "fontface.title":          "2",
    "fontface":                "1",
    "fontfamily.title":        "sans",
    "fontfamily":              "sans",
    "fontsize":                "12",
    "frame":                   "FALSE",
    "grid":                    "FALSE",
    "h":                       "-1",
    "lty.grid":                "solid",
    "lwd.border.title":        "1",
    "lwd.title":               "1",
    "lwd.grid":                "1",
    "min.distance":            "1",
    "reverseStrand":           "FALSE",
    "rotation.title":          "90",
    "showAxis":                "TRUE",
    "showTitle":               "TRUE",
    "v":                       "-1",
    "lwd.border":              "1"
}
pUcscTrack.script  = """
library (Gviz)
ucscTrack = UcscTrack (
    genome          = "{{genome}}",
    chromosome      = "{{chrom}}",
    from            = {{from}},
    to              = {{to}},
    name            = "{{name}}",
    trackType       = "{{gvizTrack}}",
    table           = "{{table}}",
    track           = "{{ucscTrack}}",
{{ proc.args | \
  (lambda x: ",\\n".join( \
     [ \
        "\\t%s\\t= %s" % (key, ('"'+ val +'"' if key in [ \
            "col.line", "col", "fill", "fontfamily.group", "fontcolor.group", "fontcolor.item", \
            "just.group", "lty", "shape", "background.panel", "background.title", "col.border.title", \
            "col.frame", "col.grid", "fontcolor.title", "fontcolor", "fontfamily.title", "fontfamily", \
            "lty.grid", "col.axis"] and val != "NULL" else val))  \
        for key, val in x.iteritems() \
     ] \
  ))(_) }}
)
saveRDS (ucscTrack, "{{outfile}}")
"""

"""
@name:
    pGenomePlot
@description:
    plot the genomic features
@input:
    `trkfiles:files`: the list of track dumped files
    `genome`:         the genome
    `chrom`:          the chromosome
    `from`:           the start position to display
    `to`:             the end position to display
@output:
    `outfile:file`:   the figure
@requires:
    [r-Gviz](https://rdrr.io/bioc/Gviz)
"""
pGenomePlot = proc ()
pGenomePlot.input  = "trkfiles:files, genome, chrom, from, to"
pGenomePlot.output = "outfile:file:genome_plot.{{#}}.png"
pGenomePlot.lang   = "Rscript"
pGenomePlot.args   = {
    "showIdeo":   True,
    "showAxis":   True,
    "resolution": 300
}
pGenomePlot.script = """
library (Gviz)
tracks = c ()
if ({{ proc.args.showIdeo | str(_).upper() }}) {
    tracks = c (tracks, IdeogramTrack(genome="{{genome}}", chromosome="{{chrom}}"))
}
if ({{ proc.args.showAxis | str(_).upper() }}) {
    tracks = c (tracks, GenomeAxisTrack())
}
for (t in noquote(unlist(strsplit("{{trkfiles | " | ".join(_)}}", " \\\\| ")))) {
    tracks = c (tracks, readRDS(t))
}
size={{proc.args.resolution}}*480/72
png (file = "{{outfile}}", width=size, height=size, res={{proc.args.resolution}})
plotTracks (
    as.list(tracks),
    from                  = {{from}},
    to                    = {{to}}
)
dev.off()
"""

