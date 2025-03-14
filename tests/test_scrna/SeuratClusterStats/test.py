import os
from biopipen.core.proc import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import (
    SeuratClustering as SeuratClustering_,
    SeuratClusterStats as SeuratClusterStats_,
    CellTypeAnnotation as CellTypeAnnotation_,
    ModuleScoreCalculator as ModuleScoreCalculator_,
    MarkersFinder,
    SeuratSubClustering as SeuratSubClustering_,
    MetaMarkers as MetaMarkers_,
    RadarPlots as RadarPlots_,
    TopExpressingGenes as TopExpressingGenes_,
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


# class PrepareSeurat(Proc):
#     """Prepare the data

#     Requires:
#         - name: Seurat
#           check: |
#             {{proc.lang}} <(echo "library(Seurat)")
#     """

#     input = "name"
#     input_data = ["pbmc3k"]
#     output = "outfile:file:{{in.name}}.RDS"
#     lang = config.lang.rscript
#     script = """
#         library(Seurat)
#         library(SeuratData)
#         set.seed(1234)
#         name <- {{in.name | r}}
#         InstallData(name)
#         data <- LoadData(name)
#         data <- UpdateSeuratObject(data)
#         data$Sample <- paste0("S", sample(1:2, ncol(data), replace = TRUE))
#         n_g1 <- floor(ncol(data) / 3)
#         data$groups <- c(rep("g1", n_g1), rep("g2", ncol(data) - n_g1))
#         data$groups <- sample(sample(data$groups))
#         data$RNA <- split(data$RNA, data$Sample)
#         # pbmc_small <- NormalizeData(pbmc_small)
#         # pbmc_small <- FindVariableFeatures(pbmc_small, selection.method = "vst", nfeatures = 2000)
#         # pbmc_small <- ScaleData(pbmc_small)
#         data <- SCTransform(data, verbose = FALSE)
#         data <- RunPCA(data, npcs = 30, verbose = FALSE)
#         # pbmc_small <- JoinLayers(pbmc_small)
#         # data <- PrepSCTFindMarkers(data)
#         saveRDS(data, {{out.outfile | quote}})
#     """  # noqa: E501


class SeuratClustering(SeuratClustering_):
    requires = PrepareSeurat
    envs = {
        "ncores": 1,
        "FindNeighbors": {"dims": 5},
        "FindClusters": {"resolution": "0.1:1,.8"},
    }


class CellTypeAnnotation(CellTypeAnnotation_):
    requires = SeuratClustering
    envs = {
        "tool": "hitype",
        "hitype_tissue": None,
        "hitype_db": "hitypedb_pbmc3k",
    }


class TopExpressingGenes(TopExpressingGenes_):
    requires = CellTypeAnnotation
    envs = {"cases": {"Cluster": {}}}


class SeuratSubClustering(SeuratSubClustering_):
    requires = CellTypeAnnotation
    envs = {
        "cases": {
            "mono_subcluster": {
                "subset": "seurat_clusters == 'FCFR3A+ Mono'",
                "FindClusters": {"resolution": "0.1:0.8,0.1"},
            },
            "dc_subcluster": {
                "subset": "seurat_clusters == 'DC'",
                "FindClusters": {"resolution": "0.1:0.8,0.2"},
            },
        }
    }


class ClusterMarkers(MarkersFinder):
    requires = SeuratSubClustering
    envs = {
        "cases": {
            # Test mixed types of cases
            "Cluster": {
                "prefix_group": False,
                "error": False,
                "allmarker_plots": {"Heatmap": {"plot_type": "heatmap"}},
            },
            "Comparison": {"group-by": "groups", "error": False, "ident-1": "g1"},
        }
    }


class DEG(MarkersFinder):
    requires = SeuratSubClustering
    envs = {
        # "mutaters": {"Cluster": "if_else(seurat_clusters %in% c('c1', 'c2', 'c3'), seurat_clusters, NA)"},
        "prefix_each": False,
        "cases": {
            "Group": {
                "group-by": "groups",
                "each": "seurat_clusters",
                # "each": "Cluster",
                "ident-1": "g1",
                "sigmarkers": "p_val < 0.5",
            }
        },
        "overlaps": {
            "Group": {
                "sigmarkers": 'abs(avg_log2FC) > 1',
                "cases": ["Group::FCFR3A+ Mono", "Group::DC", "Group::Platelet"],
                "venn": {"save_code": True},
                "upset": {"save_code": True},
            }
        },
    }
    order = 99


class MetaMarkers(MetaMarkers_):
    requires = SeuratSubClustering
    envs = {
        "group-by": "seurat_clusters",
    }


class RadarPlots(RadarPlots_):
    requires = SeuratSubClustering
    envs = {
        "by": "groups",
        "cases": {"nobreakdown": {}, "breakdown": {"breakdown": "letter.idents"}},
    }


class ModuleScoreCalculator(ModuleScoreCalculator_):
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


class SeuratClusterStats(SeuratClusterStats_):
    requires = ModuleScoreCalculator
    envs = {
        "stats": {
            "Number of cells in each cluster": {
                "pie": True,
            },
            "Number of cells in each cluster by Sample": {
                "group_by": "Sample",
                # "save_data": True,
                "frac": "group",
                "plot_type": "circos",
                "labels_rot": True,
            },
            "Number of cells in each old cluster": {
                "plot_type": "pie",
                "ident": "seurat_clusters_id",
            },
        },
        "features": {
            "Gene expressions in g1": {
                "plot_type": "ridge",
                "subset": "groups == 'g1'",
            },
            "Gene expressions in g2": {
                "plot_type": "ridge",
                "facet_ncol": 4,
                "subset": "groups == 'g2'",
            },
            "Ridge plots with ident groups": {
                "plot_type": "ridge",
                "ident": "groups",
                "features": ["CD79A", "CD79B"],
            },
            "Ridge plots with single feature": {
                # ncol automatically set to 1
                "plot_type": "ridge",
                "features": "SRSF7",
            },
            "Violin plots": {"plot_type": "violin"},
            "Violin plots (ncol=4)": {"plot_type": "violin", "facet_ncol": 4},
            "Violin plots (CD8A,NKG7)": {
                "plot_type": "violin",
                "features": ["CD8A", "NKG7"],
            },
            "Feature plot": {"plot_type": "dim", "features": "SRSF7"},
            "Dot plot": {"plot_type": "dot"},
            "Heatmap": {"plot_type": "heatmap"},
        },
        "dimplots": {
            "seurat_clusters": {"group_by": "seurat_clusters"},
            "nk_subcluster": {"group_by": "mono_subcluster"},
            "dc_subcluster": {"group_by": "dc_subcluster"},
        },
    }


def pipeline():
    # return get_pipeline(__file__, enable_report=True).set_starts(PrepareSeurat)
    return get_pipeline(__file__).set_starts(PrepareSeurat)


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
