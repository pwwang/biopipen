from pipen import Proc
from biopipen.ns.web import Download
from biopipen.ns.misc import Str2File
from biopipen.ns.gsea import FGSEA
from datar.all import flatten, tibble, select
from biopipen.core.testing import get_pipeline

FGSEA = Proc.from_proc(
    FGSEA,
    requires=[Download, Str2File],
    input_data=lambda ch1, ch2: tibble(
        *flatten(ch1),
        ch2,
        _name_repair="minimal",
    ) >> select(0, 2, 1),
    envs={"clscol": "Group", "classes": ["MMP9", "CTRL"]}
)


def pipeline():
    return get_pipeline(__file__).set_starts(
        Download,
        Str2File,
    ).set_data(
        [
            "https://www.ncbi.nlm.nih.gov/geo/download/"
            "?acc=GSE179367"
            "&format=file"
            "&file=GSE179367%5Fgene%5Fcount%2Ereal%2Etxt%2Egz",
            "https://github.com/pwwang/immunopipe-example/raw/refs/heads/master/data/"
            "KEGG_metabolism.short.gmt",
        ],
        [(
            (
                "Sample\tGroup\n"
                "shMMP9_1\tMMP9\n"
                "shMMP9_2\tMMP9\n"
                "shMMP9_3\tMMP9\n"
                "shControl_1\tCTRL\n"
                "shControl_2\tCTRL\n"
                "shControl_3\tCTRL\n"
            ),
            "samples.txt"
        )]
    )


def testing(pipen):
    # assert pipen._succeeded
    outfile = (
        pipen.procs[-1].workdir.joinpath(
            "0",
            "output",
            "acc.GSE179367.format.file.file.GSE179367_gene_count.real.txt.fgsea",
        )
    )
    assert outfile.is_dir()


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
