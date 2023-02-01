from pathlib import Path
from datar.misc import flatten

from pipen_args import args as _
from biopipen.core.proc import Proc
from biopipen.core.config import config
from biopipen.core.testing import get_pipeline
from biopipen.ns.web import Download
from biopipen.ns.scrna_metabolic_landscape import ScrnaMetabolicLandscape

HERE = Path(__file__).parent.resolve()
METADATA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/"
    "suppl/GSE123813_bcc_tcell_metadata.txt.gz"
)
COUNTDATA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/"
    "suppl/GSE123813_bcc_scRNA_counts.txt.gz"
)


class DownloadData(Download):
    """Download the data"""
    envs = {"tool": "aria2c"}

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
    envs = {
        "seed": 8525,
        "patients": ["su001", "su002"],  # None for all patients
    }
    lang = config.lang.rscript
    script = "file://scripts/PrepareData.R"


scrna_ml_pipe = ScrnaMetabolicLandscape(
    {
        "is_seurat": True,
        "gmtfile": HERE.parent.parent / "data/reference/KEGG_metabolism.gmt",
        "grouping": "seurat_clusters",
        "grouping_prefix": "Cluster",
        "subsetting": "timepoint",
        "subsetting_prefix": "Timepoint",
        "subsetting_comparison": {"post_vs_pre": ["post", "pre"]},
        "mutaters": {
            "timepoint": "treatment"
        },
    }
)

scrna_ml_pipe.procs.MetabolicExprImpute.envs["alra_args"] = {
    "use.mkl": True,
    "mkl.seed": 8525,
}


def pipeline():
    scrna_ml_pipe.procs.MetabolicInput.requires = PrepareData
    return (
        get_pipeline(__file__)
        .set_start(DownloadData)
        .set_data(
            [
                METADATA_URL,
                COUNTDATA_URL,
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
    assert outfile.is_file()


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
