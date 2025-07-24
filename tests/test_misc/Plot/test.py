from pathlib import Path
from pipen.proc import Proc
from biopipen.ns.web import Download as Download_
from biopipen.ns.misc import Plot as Plot_
from biopipen.core.config import config
from biopipen.core.testing import get_pipeline


class DownloadHeatmapData(Download_):
    """
            shControl_1     shControl_2     shControl_3     shMMP9_1        shMMP9_2        shMMP9_3
    A1BG    3.527958606     3.555128994     3.404453572     3.173882172     3.129639113     3.342048348
    A1BG-AS1        2.907038705     2.965882154     2.996968962     3.411259458     3.275175451     3.141627339
    A2M     5.199342548     5.013640911     5.182680381     6.715481584     6.743501504     6.697746024
    A2M-AS1 0.361855929     0.357423634     0.246929735     0.879261128     0.993727251     0.767122534
    """  # noqa: E501
    input_data = [
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE179367"
        "&format=file&file=GSE179367%5Fgene%5Fcount%2Ereal%2Etxt%2Egz"
    ]


class PrepareHeatmapData(Proc):
    """Convert downloaded heatmap data to a format suitable for plotting."""
    requires = DownloadHeatmapData
    input = "infile:file"
    output = "outfile:file:{{in.infile | stem}}.qs"
    lang = config.lang.rscript
    script = """
        infile <- gzfile("{{in.infile}}")
        outfile <- "{{out.outfile}}"
        data <- read.table(infile, header = TRUE, row.names = NULL, sep = "\\t")
        colnames(data)[1] <- "Gene"
        data <- data[!duplicated(data$Gene), , drop = FALSE]
        # Sample 100 rows for the heatmap
        data <- data[sample(1:nrow(data), 20), , drop = FALSE]
        biopipen.utils::save_obj(data, outfile)
    """


class HeatmapPlot(Plot_):
    requires = PrepareHeatmapData
    envs = {
        "fn": "Heatmap",
        "in_form": "wide-rows",
        "columns_by": "Gene",
        "flip": True,
    }


class PrepareManhattanData(Proc):
    """Prepare data for Manhattan plot."""
    input = "in:var"
    input_data = [1]
    output = "outfile:file:manh{{in.in}}.qs"
    lang = config.lang.rscript
    script = """
        nsim <- 50000

        simdata <- data.frame(
            "chromosome" = paste0("chr",
                                  sample(c(1:22,"X"), size = nsim, replace = TRUE)),
            "position" = sample(1:100000000, size = nsim),
            "P.value" = rbeta(nsim, shape1 = 5, shape2 = 1)^7,
            "label" = paste0("SNP", 1:nsim)
        )
        biopipen.utils::save_obj(simdata, "{{out.outfile}}")
    """


class PrepareManhattanDataRescale(Proc):
    """Prepare data for Manhattan plot with rescaling."""
    requires = PrepareManhattanData
    input = "infile:file"
    output = "outfile:file:{{in.infile | stem}}-forrescale.qs"
    lang = config.lang.rscript
    script = """
        data <- biopipen.utils::read_obj("{{in.infile}}")
        chr5_pvals <- data$P.value[data$chromosome == "chr5"]
        chr5_pvals[sample(1:length(chr5_pvals), 15)] <- exp(-rnorm(15, 50, 2))
        data$P.value[data$chromosome == "chr5"] <- chr5_pvals
        biopipen.utils::save_obj(data, "{{out.outfile}}")
    """


class ManhattanPlot(Plot_):
    """Manhattan plot using ggmanh package in R."""
    requires = PrepareManhattanData
    envs = {
        "fn": "ManhattanPlot",
        "pval_by": "P.value",
        "chr_by": "chromosome",
        "pos_by": "position",
        "title": "Simulated P.Values",
        "ylab": "P",
    }


class ManhattanPlot2(Plot_):
    """Manhattan plot with additional options."""
    requires = PrepareManhattanData
    envs = {
        "fn": "ManhattanPlot",
        "pval_by": "P.value",
        "chr_by": "chromosome",
        "pos_by": "position",
        "label_by": "label",
        "chromosomes": ["chr1", "chr2"],
        "rescale": False,
        "signif": [1e-5, 1e-3],
        "devpars": {"width": 600},
    }


class ManhattanPlotRescale(Plot_):
    """Manhattan plot with rescaling."""
    requires = PrepareManhattanDataRescale
    envs = {
        "fn": "ManhattanPlot",
        "pval_by": "P.value",
        "chr_by": "chromosome",
        "pos_by": "position",
        "label_by": "label",
        "highlight": 'chromosome == "chr5"',
        "highlight_color": "red",
        "chromosomes": ["chr1-chr10", "chrX"],
        "devpars": {"width": 600},
        "rescale": True,
        "rescale_ratio_threshold": 3,
    }


class PrepareQQPlotData(Proc):
    """Prepare data for QQ plot."""
    input = "seed:var"
    input_data = [8525]
    output = "outfile:file:data.txt"
    lang = config.lang.rscript
    script = """
        set.seed({{in.seed}})
        data <- data.frame(norm = rnorm(100))
        write.table(data,
                    file = "{{out.outfile}}",
                    sep = "\\t",
                    quote = FALSE,
                    row.names = FALSE)
    """


class QQPlot(Plot_):
    """QQ plot using plotthis package in R."""
    requires = PrepareQQPlotData
    envs = {
        "fn": "QQPlot",
        "val": "norm",
        "band": True,
        "read_opts": {"header": True, "sep": "\t"},
    }


class ROCPlotSingle(Plot_):
    """ROC plot using plotthis package in R."""
    input_data = [Path(__file__).parent.joinpath("data", "single.txt")]
    envs = {
        "fn": "ROCCurve",
        "truth_by": "D",
        "score_by": "M1",
        "read_opts": {"header": True, "sep": "\t"},
    }


class ROCPlotMulti(Plot_):
    """ROC plot with multiple scores."""
    input_data = [Path(__file__).parent.joinpath("data", "multi.txt")]
    envs = {
        "fn": "ROCCurve",
        "truth_by": "D",
        "score_by": ["M1", "M2"],
        "read_opts": {"header": True, "sep": "\t"},
    }


def pipeline():
    return get_pipeline(__file__).set_starts(
        DownloadHeatmapData,
        PrepareManhattanData,
        PrepareQQPlotData,
        ROCPlotSingle,
        ROCPlotMulti,
    )


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
