from pathlib import Path
from pipen import Proc  # type: ignore
from biopipen.core.config import config
from biopipen.ns.scrna import (
    SeuratTo10X as SeuratTo10X_,
    SeuratPreparing as SeuratPreparing_,
    PseudoBulkDEG as PseudoBulkDEG_,
    SeuratClustering as SeuratClustering_,
    CellTypeAnnotation as CellTypeAnnotation_,
    TopExpressingGenes as TopExpressingGenes_,
    SeuratSubClustering as SeuratSubClustering_,
    MarkersFinder,
    ModuleScoreCalculator as ModuleScoreCalculator_,
    SeuratClusterStats as SeuratClusterStats_,
    ScFGSEA as ScFGSEA_,
)
from biopipen.core.testing import get_pipeline


class PrepareSeurat(Proc):
    """Prepare the data"""

    input = "name"
    input_data = ["pbmc3k"]
    output = "outfile:file:{{in.name}}.RDS"
    lang = config.lang.rscript
    script = """
        set.seed(8525)
        options(timeout=600)
        library(SeuratData)
        InstallData("pbmc3k")
        pbmc3k <- Seurat::UpdateSeuratObject(pbmc3k)
        pbmc3k$Sample <- paste0("S", sample(1:12, nrow(pbmc3k), replace = TRUE))
        saveRDS(pbmc3k, {{out.outfile | quote}})
    """


class SeuratTo10X(SeuratTo10X_):
    requires = PrepareSeurat
    envs = {"split_by": "Sample"}


class PrepareMetafile(Proc):
    """Prepare the metafile for SeuratPreparing"""

    requires = SeuratTo10X
    input = "indir:dir"
    output = "metafile:file:SeuratPreparing-test.txt"
    lang = config.lang.python
    script = """
        import os
        indir = {{in.indir | quote}}
        metafile = {{out.metafile | quote}}
        with open(metafile, "w") as f:
            f.write("Sample\\tGroup\\tEach\\tRNAData\\n")
            for s in os.listdir(indir):
                if os.path.isdir(os.path.join(indir, s)):
                    group = (
                        "Group1"
                        if s in ["S1", "S2", "S3", "S4", "S5", "S6"]
                        else "Group2"
                    )
                    each = (
                        "E1"
                        if s in ["S1", "S2", "S3", "S7", "S8", "S9"]
                        else "E2"
                    )
                    f.write(f"{s}\\t{group}\\t{each}\\t{os.path.join(indir, s)}\\n")
    """


class SeuratPreparing(SeuratPreparing_):
    requires = PrepareMetafile
    envs = {
        "cell_qc": "runif(n()) < 0.5",
        "IntegrateLayers": {"method": "rpca", "k-weight": 30},
    }


class SeuratPreparing2(SeuratPreparing_):
    requires = PrepareMetafile
    envs = {
        "cell_qc": "runif(n()) < 0.5",
        "cell_qc_per_sample": True,
        "doublet_detector": "DoubletFinder",
        "DoubletFinder": {"PCs": 3},
    }


class SeuratPreparing3(SeuratPreparing_):
    requires = PrepareMetafile
    envs = {
        "mutaters": {"X": "1"},
        "cell_qc": "runif(n()) < 0.5",
        "cell_qc_per_sample": True,
        "doublet_detector": "scDblFinder",
    }


class SeuratClustering(SeuratClustering_):
    requires = SeuratPreparing
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


class CellTypeAnnotationDirect(CellTypeAnnotation_):
    requires = SeuratClustering
    envs = {
        "tool": "direct",
        "cell_types": [
            "Naive CD4+ T",
            "B",
            "Memory CD4+",
            "Naive CD4+ T",
            "DC",
            "DC",
            "CD8+ T",
            "NK",
            "FCFR3A+ Mono",
            "CD8+ T",
        ],
        "merge": True,
    }


class TopExpressingGenes(TopExpressingGenes_):
    requires = CellTypeAnnotation
    envs = {"cases": {"Cluster": {}}}


class SeuratSubClustering(SeuratSubClustering_):
    requires = SeuratClustering
    envs = {
        "cache": False,
        "cases": {
            "mono_subcluster": {
                "subset": "seurat_clusters == 'c1'",
                "FindClusters": {"resolution": "0.1:0.8,0.1"},
            },
            "dc_subcluster": {
                "subset": "seurat_clusters == 'c2'",
                "FindClusters": {"resolution": "0.1:0.8,0.2"},
            },
        },
    }


class ClusterMarkers(MarkersFinder):
    requires = SeuratSubClustering
    envs = {
        "cases": {
            # Test mixed types of cases
            "Cluster": {
                "error": False,
                "sigmarkers": "p_val_adj < 0.05 & avg_log2FC > 0",
                "marker_plots": {
                    "Heatmap": {"plot_type": "heatmap"},
                },
                "allmarker_plots": {
                    "Heatmap": {
                        "plot_type": "heatmap",
                    },
                },
            },
            "Comparison": {
                "sigmarkers": "p_val < 0.1",
                "group_by": "Group",
                "error": False,
                "ident_1": "Group1",
            },
        }
    }


class DEGSingleComparison(MarkersFinder):
    requires = SeuratSubClustering
    envs = {
        "group_by": "Group",
        "ident_1": "Group1",
        "ident_2": "Group2",
        "sigmarkers": "p_val < 0.5",
    }


class DEGSingleComparisonWithEach(MarkersFinder):
    requires = SeuratSubClustering
    envs = {
        "group_by": "Group",
        "each": "seurat_clusters",
        "ident_1": "Group1",
        "sigmarkers": "p_val < 0.5",
        "allmarker_plots": {
            "Heatmap": {
                "plot_type": "heatmap",
            },
            "Dot": {
                "plot_type": "dot",
            },
        },
    }


class DEG(MarkersFinder):
    requires = SeuratSubClustering
    envs = {
        # "mutaters": {"Cluster": "if_else(seurat_clusters %in% c('c1', 'c2', 'c3'), seurat_clusters, NA)"},
        "cases": {
            "Group": {
                "group_by": "Group",
                "each": "seurat_clusters.0.8",
                "subset": "seurat_clusters.0.8 %in% c('c1', 'c2', 'c3')",
                # "each": "Cluster",
                "ident_1": "Group1",
                "sigmarkers": "p_val < 0.5",
                "overlaps": {
                    "Venn": {
                        "sigmarkers": "abs(avg_log2FC) > 1",
                        "plot_type": "venn",
                    },
                },
            },
        },
    }
    order = 99


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
                "ident": "seurat_clusters.0.2",
            },
        },
        "features": {
            "Gene expressions in g1": {
                "plot_type": "ridge",
                "subset": "Group == 'Group1'",
            },
            "Gene expressions in g2": {
                "plot_type": "ridge",
                "facet_ncol": 4,
                "subset": "Group == 'Group2'",
            },
            "Ridge plots with ident groups": {
                "plot_type": "ridge",
                "ident": "Group",
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
            "nk_subcluster": {
                "group_by": "mono_subcluster",
                "reduction": "mono_subcluster.umap",
            },
            "dc_subcluster": {
                "group_by": "dc_subcluster",
                "reduction": "dc_subcluster.umap",
            },
        },
    }


class PseudoBulkDEG(PseudoBulkDEG_):
    requires = SeuratPreparing
    envs = {
        "group_by": "Group",
        "ident_1": "Group1",
        "ident_2": "Group2",
        "ncores": 2,
        "sigmarkers": "p_val < 0.25",
        "error": False,
        "aggregate_by": ["Sample", "Group"],
    }


class PseudoBulkDEGEach(PseudoBulkDEG_):
    requires = SeuratPreparing
    envs = {
        "group_by": "Group",
        "ident_1": "Group1",
        "ident_2": "Group2",
        "each": "Each",
        "sigmarkers": "p_val < 0.25",
        "error": False,
        "aggregate_by": ["Sample", "Group"],
        "allmarker_plots": {
            "Heatmap": {
                "plot_type": "heatmap",
            },
        },
        "allenrich_plots": {
            "Heatmap": {
                "plot_type": "heatmap",
            },
        },
        "overlaps": {
            "Venn": {
                "sigmarkers": "abs(log2FC) > 1",
                "plot_type": "venn",
            },
        },
    }


class CellTypeAnnotationSCCatch(CellTypeAnnotation_):
    requires = SeuratClustering
    envs = {
        "tool": "sccatch",
        "sccatch_args": {
            "species": "Human",
            "tissue": "Peripheral blood",
        },
    }


class CellTypeAnnotationSCCatch2(CellTypeAnnotation_):
    requires = SeuratClustering
    envs = {
        "tool": "sccatch",
        "merge": True,
        "sccatch_args": {
            "marker": str(Path(__file__).parent.joinpath("data", "tcell.sccatch.RDS")),
        },
    }


class SeuratClusterStatsSCCatch(SeuratClusterStats_):
    requires = CellTypeAnnotationSCCatch
    envs = {
        "stats": {
            "Number of cells in each cluster by Sample": {
                "group_by": "seurat_clusters",
            }
        }
    }


class SeuratClusterStatsSCCatch2(SeuratClusterStats_):
    requires = CellTypeAnnotationSCCatch2
    envs = {
        "stats": {
            "Number of cells in each cluster by Sample": {
                "group_by": "seurat_clusters",
            }
        }
    }


class ScFGSEASingle(ScFGSEA_):
    requires = SeuratClustering
    envs = {
        "ident_1": "Group2",
        "ident_2": "Group1",
        "group_by": "Group",
        "subset": "seurat_clusters == 'c1'",
        "minSize": 1,
        "gmtfile": (
            "https://raw.githubusercontent.com/pwwang/immunopipe-example/"
            "master/data/KEGG_metabolism.short.gmt"
        ),
    }


class ScFGSEAEach(ScFGSEA_):
    requires = SeuratClustering
    envs = {
        # Test if seurat_clusters == 1 gets ignored
        "mutaters": {"Group2": "ifelse(seurat_clusters == 'c1', NA, Group)"},
        "ident_1": "Group1",
        "ident_2": "Group2",
        "group_by": "Group2",
        "minsize": 1,
        # "gmtfile": Path(__file__).parent.parent.parent.joinpath(
        #     "data/reference/KEGG_metabolism.gmt"
        # ),
        "each": "seurat_clusters",
        "gmtfile": (
            "https://raw.githubusercontent.com/pwwang/immunopipe-example/"
            "master/data/KEGG_metabolism.short.gmt"
        ),
    }


def pipeline():
    return (
        get_pipeline(__file__)
        # get_pipeline(__file__, enable_report=True)
        .set_starts(PrepareSeurat)
    )


def testing(pipen):
    # assert pipen._succeeded
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
