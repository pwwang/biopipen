
from biopipen.core.config import config
from biopipen.core.proc import Proc
from biopipen.ns.plot import Manhattan as Manhattan_
from biopipen.core.testing import get_pipeline


class PrepareData(Proc):
    """Prepare the data for plotting"""
    input = "in:var"
    output = "outfile:file:manh{{in.in}}.txt"
    lang = config.lang.rscript
    script = """
        library(ggmanh)
        library(SeqArray)

        set.seed(1000)
        nsim <- 50000

        simdata <- data.frame(
            "chromosome" = paste0("chr",
                                  sample(c(1:22,"X"), size = nsim, replace = TRUE)),
            "position" = sample(1:100000000, size = nsim),
            "P.value" = rbeta(nsim, shape1 = 5, shape2 = 1)^7,
            "label" = paste0("SNP", 1:nsim)
        )

        write.table(simdata,
                    file = "{{out.outfile}}",
                    sep = "\\t",
                    quote = FALSE,
                    row.names = FALSE)
    """


class Manhattan(Manhattan_):
    requires = PrepareData
    envs = {"rescale": False}


class Manhattan2(Manhattan_):
    requires = PrepareData
    envs = {
        "label_col": "label",
        "zoom": ["chr1", "chr2"],
        "rescale": False,
        "signif": [1e-5, 1e-3],
    }


class PrepareForRescale(Proc):
    requires = PrepareData
    input = "infile:file"
    output = "outfile:file:{{in.infile | stem}}-forrescale.txt"
    lang = config.lang.rscript
    script = """
        data <- read.table(
            "{{in.infile}}", header = TRUE, sep = "\\t", check.names = FALSE)
        chr5_pvals <- data$P.value[data$chromosome == "chr5"]
        chr5_pvals[sample(1:length(chr5_pvals), 15)] <- exp(-rnorm(15, 50, 2))
        data$P.value[data$chromosome == "chr5"] <- chr5_pvals
        write.table(
            data,
            file = "{{out.outfile}}", sep = "\\t", quote = FALSE, row.names = FALSE
        )
    """


class Manhattan3(Manhattan_):
    requires = PrepareForRescale
    envs = {
        "label_col": "label",
        "hicolors": "red",
        "chroms": ["chr1-chr10", "chrX"],
        "devpars": {"width": 600},
        "rescale": True,
        "rescale_ratio_threshold": 3,
    }


def pipeline():
    return get_pipeline(__file__).set_starts(
        PrepareData
    ).set_data(
        [1]
    )


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
