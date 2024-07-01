
from biopipen.core.config import config
from biopipen.core.proc import Proc
from biopipen.ns.plot import QQPlot as QQPlot_
from biopipen.core.testing import get_pipeline


class PrepareData(Proc):
    """Prepare the data for plotting"""
    input = "seed:var"
    output = "outfile:file:data.txt, theorfile:file:theor.txt"
    lang = config.lang.rscript
    script = """
        set.seed({{in.seed}})
        smp <- data.frame(norm = rnorm(100))
        write.table(smp,
                    file = "{{out.outfile}}",
                    sep = "\\t",
                    quote = FALSE,
                    row.names = FALSE)
        theor <- data.frame(theor = exp(rnorm(200)))
        write.table(theor,
                    file = "{{out.theorfile}}",
                    sep = "\\t",
                    quote = FALSE,
                    row.names = FALSE)
    """


class QQPlotDefault(QQPlot_):
    requires = PrepareData


class QQPlotCustom(QQPlot_):
    requires = PrepareData
    envs = {
        "theor_col": "theor",
        "theor_trans": "log",
        "args": {"distribution": "custom"},
        "band": {"disabled": True},
    }


def pipeline():
    return get_pipeline(__file__).set_starts(PrepareData).set_data([0])


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
