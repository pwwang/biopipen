from biopipen.core.proc import Proc
from biopipen.core.config import config
from biopipen.core.testing import get_pipeline
from biopipen.ns.web import Download
from biopipen.ns.scrna_metabolic_landscape import ScrnaMetabolicLandscape

METADATA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl/"
    "GSE123813%5Fbcc%5Ftcell%5Fmetadata.txt.gz"
)
COUNTDATA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl/"
    "GSE123813%5Fbcc%5FscRNA%5Fcounts.txt.gz"
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
    input_data = lambda ch: [tuple(ch.iloc[:, 0])]
    output = "outfile:file:{{in.metafile | stem}}.RDS"
    envs = {
        "seed": 8525,
        "patients": ["su001", "su002"],  # None for all patients
    }
    lang = config.lang.rscript
    script = "file://scripts/PrepareData.R"


scrna_ml_pipe = ScrnaMetabolicLandscape(
    is_seurat=True,
    gmtfile=(
        "https://raw.githubusercontent.com/pwwang/immunopipe-example/master/data/"
        "KEGG_metabolism.short.gmt"
    ),
    group_by="cluster",
    mutaters={"timepoint": "treatment"},
)

if "MetabolicExprImputation" in scrna_ml_pipe.procs:
    scrna_ml_pipe.procs.MetabolicExprImputation.envs["alra_args"] = {
        "use.mkl": False,
        "mkl.seed": 8525,
    }


scrna_ml_pipe.procs.MetabolicPathwayActivity.envs["cases"] = {
    "NoSubsetting": {},
    "ByTimepoint": {
        "subset_by": "timepoint",
        "plots": {
            "Merged Heatmap": {"plot_type": "merged_heatmap"},
        },
    },
}


def pipeline():
    scrna_ml_pipe.procs.MetabolicInput.requires = PrepareData
    return (
        get_pipeline(__file__)
        # get_pipeline(__file__, enable_report=True)
        .set_start(DownloadData)
        .set_data(
            [
                METADATA_URL,
                COUNTDATA_URL,
            ]
        )
    )


def testing(pipen):
    # assert pipen._succeeded
    # outfile = pipen.outdir.joinpath(
    #     "REPORTS",
    #     "index.html",
    # )
    # assert outfile.is_file()
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
