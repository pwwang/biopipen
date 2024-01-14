from pipen import Proc
from biopipen.core.config import config
from biopipen.ns.tcr import Immunarch
from biopipen.core.testing import get_pipeline


class PrepareImmdata(Proc):
    """Prepare the data
    """

    input = "name"
    input_data = ["immdata"]
    output = "outfile:file:{{in.name}}.RDS"
    lang = config.lang.rscript
    script = """
        library(immunarch)
        data(immdata)
        immdata$data <- lapply(immdata$data, function(x) {
            x$Barcode <- sapply(seq_along(x$CDR3.aa), function(i) {
                paste(x$CDR3.aa[i], 1:x$Clones[i], sep = "_", collapse = ";")
            })
            x
        })
        saveRDS(immdata, {{out.outfile | quote}})
    """


class Immunarch(Immunarch):
    requires = PrepareImmdata
    envs = {
        "divs": {
            "by": "Status,Sex",
            "cases": {
                "chao1_bar": {"method": "chao1", "plot_type": "bar"},
                "chao1_box": {"method": "chao1", "plot_type": "box"},
                "hill_bar": {"method": "hill", "plot_type": "bar"},
                "hill_box": {"method": "hill", "plot_type": "box"},
                "div_bar": {"method": "div", "plot_type": "bar"},
                "div_box": {"method": "div", "plot_type": "box"},
                "gini.simp_bar": {"method": "gini.simp", "plot_type": "bar"},
                "gini.simp_box": {"method": "gini.simp", "plot_type": "box"},
                "inv.simp_bar": {"method": "inv.simp", "plot_type": "bar"},
                "inv.simp_box": {"method": "inv.simp", "plot_type": "box"},
                "gini_bar": {"method": "gini", "plot_type": "bar"},
                "gini_box": {"method": "gini", "plot_type": "box"},
                "d50_bar": {"method": "d50", "plot_type": "bar"},
                "d50_box": {"method": "d50", "plot_type": "box"},
            },
        }
    }


def pipeline():
    return get_pipeline(__file__, plugins=["no:report"]).set_starts(
        PrepareImmdata
    )


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
