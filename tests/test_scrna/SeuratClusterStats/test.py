from pipen import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import (
    SeuratClustering,
    SeuratClusterStats,
    CellTypeAnnotation,
    ModuleScoreCalculator,
    MarkersFinder,
    SeuratSubClustering,
    MetaMarkers,
    RadarPlots,
    TopExpressingGenes,
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
        pbmc_small$Sample <- pbmc_small$letter.idents
        pbmc_small$RNA <- split(pbmc_small$RNA, pbmc_small$Sample)
        # pbmc_small <- NormalizeData(pbmc_small)
        # pbmc_small <- FindVariableFeatures(pbmc_small, selection.method = "vst", nfeatures = 2000)
        # pbmc_small <- ScaleData(pbmc_small)
        pbmc_small <- SCTransform(pbmc_small, verbose = FALSE)
        pbmc_small <- RunPCA(pbmc_small, npcs = 30, verbose = FALSE)
        # pbmc_small <- JoinLayers(pbmc_small)
        pbmc_small <- PrepSCTFindMarkers(pbmc_small)
        saveRDS(pbmc_small, {{out.outfile | quote}})
    """  # noqa: E501


class SeuratClustering(SeuratClustering):
    requires = PrepareSeurat
    envs = {
        "ncores": 1,
        "FindNeighbors": {"dims": 5},
        "FindClusters": {"resolution": "0.5, 0.8"},
    }


class CellTypeAnnotation(CellTypeAnnotation):
    requires = SeuratClustering
    envs = {
        "tool": "hitype",
        "hitype_tissue": None,
        "hitype_db": "hitypedb_pbmc3k",
    }


class TopExpressingGenes(TopExpressingGenes):
    requires = CellTypeAnnotation
    envs = {"cases": {"Cluster": {}}}


class SeuratSubClustering(SeuratSubClustering):
    requires = CellTypeAnnotation
    envs = {
        "cases": {
            "mono_subcluster": {
                "subset": "seurat_clusters == 'FCFR3A+ Mono'",
                "FindClusters": {"resolution": "0.5,0.8"},
            },
            "dc_subcluster": {
                "subset": "seurat_clusters == 'DC'",
                "FindClusters": {"resolution": "0.5,0.8"},
            },
        }
    }


class ClusterMarkers(MarkersFinder):
    requires = SeuratSubClustering
    envs = {
        "cases": {
            # Test mixed types of cases
            "Cluster": {"prefix_group": False},
            "Comparison": {"group-by": "groups", "ident-1": "g1"},
        }
    }


class DEG(MarkersFinder):
    requires = SeuratSubClustering
    envs = {
        "prefix_each": False,
        "cases": {
            "Group": {
                "group-by": "groups",
                "each": "seurat_clusters",
                "ident-1": "g1",
            }
        },
        "overlap": {"Group": {}},
    }
    order = 99


class MetaMarkers(MetaMarkers):
    requires = SeuratSubClustering
    envs = {
        "group-by": "seurat_clusters",
    }


class RadarPlots(RadarPlots):
    requires = SeuratSubClustering
    envs = {"by": "groups"}


class ModuleScoreCalculator(ModuleScoreCalculator):
    requires = SeuratSubClustering
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
                "circos": True,
                "circos_labels_rot": True,
            },
            "Number of cells in each old cluster": {
                "pie": True,
                "ident": "seurat_clusters_id"
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
        },
        "dimplots": {
            "seurat_clusters": {},
            "nk_subcluster": {"ident": "mono_subcluster"},
            "dc_subcluster": {"ident": "dc_subcluster"},
        }
    }


def pipeline():
    return get_pipeline(__file__).set_starts(PrepareSeurat)


def testing(pipen):
    assert pipen._succeeded
    outfile = (
        pipen.procs[-1].workdir.joinpath(
            "0",
            "output",
            "pbmc_small.markers/OVERLAPPING/Group/markers.txt",
        )
    )
    assert outfile.is_file(), str(outfile)


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
