from pipen import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import (
    SeuratClustering,
    SeuratClusterStats,
    CellTypeAnnotation,
    ModuleScoreCalculator,
    MarkersFinder,
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
    envs = {"ncores": 1}


class CellTypeAnnotation(CellTypeAnnotation):
    requires = SeuratClustering
    envs = {
        "tool": "hitype",
        "hitype_tissue": None,
        "hitype_db": "hitypedb_pbmc3k",
    }


class MarkersFinder(MarkersFinder):
    requires = CellTypeAnnotation


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
        "stats": {
            "Number of cells in each cluster": {
                "pie": True,
            },
            "Number of cells in each cluster by Sample": {
                "group-by": "Sample",
                "table": True,
                "frac": True,
            },
            "Number of cells in each old cluster": {
                "pie": True,
                "ident": "seurat_clusters_id",
            },
        },
        "features": {
            "Gene expressions in g1": {
                "kind": "ridge",
                "subset": "groups == 'g1'",
            },
            "Gene expressions in g2": {
                "kind": "ridge",
                "ncol": 4,
                "subset": "groups == 'g2'",
            },
            "Ridge plots with ident groups": {
                "kind": "ridge",
                "ident": "groups",
                "features": "CD1C,RGS1",
            },
            "Ridge plots with single feature": {
                # ncol automatically set to 1
                "kind": "ridge",
                "features": "SRSF7",
                "plus": "theme_gray(base_size=10)",
            },
            "Violin plots": {"kind": "violin", "pt-size": 0},
            "Feature plot": {"kind": "feature", "features": "SRSF7"},
            "Dot plot": {"kind": "dot", "plus": "RotatedAxis()"},
            "Heatmap": {"kind": "heatmap"},
            "Gene expression table": {"kind": "table"},
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
            "pbmc_small.markers/DEFAULT/B/markers.txt",
        )
    )
    assert outfile.is_file(), str(outfile)


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
