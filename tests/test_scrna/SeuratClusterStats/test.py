from pipen import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import (
    SeuratClustering,
    SeuratClusterStats,
    CellTypeAnnotation,
    ModuleScoreCalculator,
)
from biopipen.core.testing import get_pipeline


class PrepareSeurat(Proc):
    """Prepare the data

    Requires:
        - name: Seurat
          check: |
            {{proc.lang}} <(echo "library(Seurat)")
    """

    input = "name"
    input_data = ["pbmc_small"]
    output = "outfile:file:{{in.name}}.RDS"
    lang = config.lang.rscript
    script = """
        library(Seurat)
        pbmc_small$Sample = pbmc_small$letter.idents
        saveRDS(pbmc_small, {{out.outfile | quote}})
    """


class SeuratClustering(SeuratClustering):
    requires = PrepareSeurat
    envs = {
        "FindIntegrationAnchors": {"reduction": "cca"},
        "IntegrateData": {"k-weight": 5},
    }


class CellTypeAnnotation(CellTypeAnnotation):
    requires = SeuratClustering
    envs = {
        "tool": "hitype",
        "hitype_tissue": None,
        "hitype_db": "hitypedb_pbmc3k",
    }


class ModuleScoreCalculator(ModuleScoreCalculator):
    requires = CellTypeAnnotation
    envs = {
        "modules": {
            # "CellCycle": {"features": "cc.genes.updated.2019"},
            # "Exhaustion": {"features": "HAVCR2,ENTPD1,LAYN,LAG3"},
            # "Activation": {"features": "IFNG"},
            # "Proliferation": {"features": "STMN1,TUBB"},
            "SomeModule": {"features": "CD3D,GZMM,CD8A,GNLY", "ctrl": 4},
        }
    }


class SeuratClusterStats(SeuratClusterStats):
    requires = ModuleScoreCalculator
    envs = {
        "features": {
            "ridgeplots_1": {
                "title": "Gene expressions in g1",
                "subset": "groups == 'g1'",
            },
            "ridgeplots_2": {
                "title": "Gene expressions in g2",
                "subset": "groups == 'g2'",
            },
            "vlnplots": {"boxplot": {}, "pt-size": 0},
            "dotplot": {"plus": "RotatedAxis()"},
            "heatmap": {"downsample": "average"},
        }
    }


def pipeline():
    # return get_pipeline(__file__).set_starts(
    return get_pipeline(__file__, plugins=["no:report"]).set_starts(
        PrepareSeurat
    )


def testing(pipen):
    outfile = (
        pipen.procs[-1].workdir.joinpath(
            "0",
            "output",
            "pbmc_small.annotated.cluster_stats",
            "features",
            "ridgeplots-1.png",
        )
    )
    assert outfile.is_file(), str(outfile)


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
