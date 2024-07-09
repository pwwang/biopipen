
from biopipen.core.config import config
from biopipen.core.proc import Proc
from biopipen.ns.plot import Scatter as Scatter_
from biopipen.core.testing import get_pipeline


class PrepareData(Proc):
    """Prepare the data for plotting"""
    input = "seed:var"
    output = "outfile:file:data.txt"
    lang = config.lang.rscript
    script = """
        set.seed({{in.seed}})
        x <- (1:100) / 10
        y <- x + rnorm(length(x))
        my.data <- data.frame(X = x,
                              Y = y,
                              y.desc = - y,
                              group = c("A", "B"))
        write.table(my.data,
                    file = "{{out.outfile}}",
                    sep = "\\t",
                    quote = FALSE,
                    row.names = FALSE)
    """


class Scatter(Scatter_):
    requires = PrepareData
    envs = {"formula": "y ~ x"}


class ScatterGroup(Scatter_):
    requires = PrepareData
    envs = {
        "formula": "y ~ x",
        "mapping": "color = group",
        "stats": {
            "poly_line": {},
            "poly_eq": {
                "mapping": """
                    aes(
                        label = paste(
                            after_stat(eq.label),
                            after_stat(rr.label),
                            sep = '*", "*'
                        )
                    )
                """
            },
        }
    }


class ScatterFacet(Scatter_):
    requires = PrepareData
    envs = {
        "formula": "y ~ x",
        "ggs": "facet_wrap(~ group)",
        "devpars": {"width": 1400},
        "stats": {
            "poly_line": {},  # formula added automatically
            "poly_eq": {"formula": "formula"},
        }
    }


def pipeline():
    return get_pipeline(__file__).set_starts(PrepareData).set_data([4321])


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
