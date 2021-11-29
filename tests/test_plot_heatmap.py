from pipen import Pipen, Proc
from pipen_args import args
from biopipen.namespaces.web import Download
from biopipen.namespaces.plot import Heatmap
from datar.all import tibble, flatten

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

pipen = (
    Pipen("TestPlotHeatmap")
    .set_start(Download)
    .set_data(
        [
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE179367&format=file&file=GSE179367%5Fgene%5Fcount%2Ereal%2Etxt%2Egz",
        ]
    )
)

if __name__ == "__main__":
    pipen.run()
