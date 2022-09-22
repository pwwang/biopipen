
from datar.all import tibble, flatten
from pipen import Proc
from biopipen.ns.web import Download
from biopipen.ns.plot import Heatmap
from biopipen.core.testing import get_pipeline


Heatmap = Proc.from_proc(
    Heatmap,
    requires=Download,
    input_data=lambda ch: tibble(infile=ch, annofiles=[flatten(ch)]),
    envs={
        "globals": "data = head(data, 100)",
        "args": {
            "right_annotation": """r:rowAnnotation(
                Boxplot = anno_boxplot(as.matrix(head(annos, 100)), outline = F)
            )"""
        }
    }
)


def pipeline():
    return get_pipeline(__file__, plugins=["no:report"]).set_starts(
        Download
    ).set_data(
        [
            "https://www.ncbi.nlm.nih.gov/geo/download/"
            "?acc=GSE179367"
            "&format=file"
            "&file=GSE179367%5Fgene%5Fcount%2Ereal%2Etxt%2Egz",
        ]
    )


def testing(pipen):
    outfile = (
        pipen.procs[-1].workdir.joinpath(
            "0",
            "output",
            "acc.heatmap/acc.heatmap.png",
        )
    )
    assert outfile.is_file()


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
