from pathlib import Path
from datar.tibble import tibble
from datar.datar import flatten

from pipen import Proc
from pipen_args import args as _
from biopipen.core.config import config
from biopipen.core.testing import get_pipeline
from biopipen.ns.web import Download
from biopipen.ns.scrna_metabolic import build_processes


class DownloadData(Download):
    """Download the data"""


class PrepareData(Proc):
    """Prepare the data

    Requires:
        - name: r-seurat
          check: |
            {{proc.lang}} <(echo "library(Seurat)")
    """
    requires = DownloadData
    input = "metafile:file, countfile:file"
    input_data = lambda ch: [tuple(flatten(ch))]
    output = "outfile:file:{{in.metafile | stem}}.RDS"
    envs = {"seed": 8525}
    lang = config.lang.rscript
    script = "file://scripts/PrepareData.R"


MetabolicInputs = build_processes({"clustered": True})
MetabolicInputs.requires = PrepareData
MetabolicInputs.input_data = lambda ch: tibble(
    metafile=ch,
    gmtfile=Path(__file__).parent.parent.parent.joinpath(
        "data/KEGG_metabolism.gmt"
    ),
    config=[
        """
        [grouping]
        groupby = "seurat_clusters"

        [subsetting]
        alias = "Timepoint"
        groupby = "timepoint"

        [subsetting.mutaters]
        timepoint = "if_else(patient != 'su001', NA_character_, treatment)"

        [design]
        post_vs_pre = ["post", "pre"]
        """
    ],
)


def pipeline():
    return (
        get_pipeline(__file__)
        .set_start(DownloadData)
        .set_data(
            [
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/"
                "suppl/GSE123813_bcc_tcell_metadata.txt.gz",
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/"
                "suppl/GSE123813_bcc_scRNA_counts.txt.gz",
            ]
        )
    )


def testing(pipen):
    outfile = (
        pipen.outdir.joinpath(
            "REPORTS",
            "index.html",
        )
    )
    assert outfile.is_dir()


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
