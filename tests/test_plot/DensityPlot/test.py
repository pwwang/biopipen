from biopipen.core.config import config
from biopipen.core.proc import Proc
from biopipen.ns.plot import DensityPlot as DensityPlot_
from biopipen.core.testing import get_pipeline


class PrepareData(Proc):
    input = "seed:var"
    output = "outfile:file:{{in.seed}}.density_data.txt"
    lang = config.lang.rscript
    script = """
        seed <- {{in.seed | r}}
        outfile <- {{out.outfile | r}}
        data <- data.frame(
            x = c(rnorm(500, -1), rnorm(500, 1)),
            group = rep(c("A", "B"), each = 500),
            facet = sample(c("F1", "F2"), 1000, replace = TRUE)
        )
        write.table(data, file=outfile, sep="\\t", row.names=FALSE, quote=FALSE)
    """


class DensityPlot1(DensityPlot_):
    requires = [PrepareData]
    envs = {
        "val_col": "x",
    }


class DensityPlot2(DensityPlot_):
    requires = [PrepareData]
    envs = {
        "val_col": "x",
        "group_by": "group",
        "facet_by": "facet",
    }


def pipeline():
    return (
        get_pipeline(__file__)
        .set_starts(PrepareData)
        .set_data([8525])
    )


def testing(pipen):
    # assert pipen._succeeded
    outfile = (
        pipen.procs[-1].workdir.joinpath("0", "output", "8525.density_data.density.png")
    )
    assert outfile.is_file()


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
