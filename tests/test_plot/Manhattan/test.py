
from datar.all import tibble, flatten
from biopipen.core.config import config
from biopipen.core.proc import Proc
from biopipen.ns.web import Download
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


class Manhattan2(Manhattan_):
    requires = PrepareData
    envs = {
        "label_col": "label",
        "zoom": ["chr1", "chr2"],
        "signif": [1e-5, 1e-3],
    }


class Manhattan3(Manhattan_):
    requires = PrepareData
    envs = {
        "label_col": "label",
        "hicolors": "red",
        "signif": [1e-5, 1e-3],
        "chroms": ["chr1-chr10", "chrX"],
        "devpars": {"width": 600},
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
