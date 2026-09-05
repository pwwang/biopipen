"""Tests for the HdWGCNA process, replicating the hdWGCNA tutorials.

Local-only test (see run.env): downloads the Zhou et al. 2020 control
snRNA-seq dataset (~1GB) and runs full WGCNA network construction.
Run with: python tests/test_scrna/HdWGCNA/test.py

Each test class replicates one hdWGCNA tutorial:
- TestHdWGCNAMetacells:  basic tutorial (harmony metacells, INH network)
  https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html
- TestHdWGCNAPseudobulk: pseudobulk tutorial (VST pseudobulk network)
  https://smorabit.github.io/hdWGCNA/articles/pseudobulk.html
- TestHdWGCNAConsensus:  consensus network (pseudobulk, split by msex)
  https://smorabit.github.io/hdWGCNA/articles/consensus_wgcna.html
- TestHdWGCNADownstream: DME analysis + module-trait correlation
  https://smorabit.github.io/hdWGCNA/articles/differential_MEs.html
  https://smorabit.github.io/hdWGCNA/articles/module_trait_correlation.html
- TestHdWGCNAEnrichment: module enrichment
  https://smorabit.github.io/hdWGCNA/articles/enrichment_analysis.html
- TestHdWGCNAViz:        network visualizations
  https://smorabit.github.io/hdWGCNA/articles/network_visualizations.html
- TestHdWGCNAPreservation: module preservation (query == ref)
  https://smorabit.github.io/hdWGCNA/articles/module_preservation.html
- TestHdWGCNAMotifs:     motif overlap analysis
  https://smorabit.github.io/hdWGCNA/articles/motif_analysis.html
- TestHdWGCNA_TFNetwork: TF regulatory network
  https://smorabit.github.io/hdWGCNA/articles/tf_network.html

Not covered (documented): spatial (needs ST data), other_metacells (needs
custom metacell reductions), projection (needs snATAC data), PPI (needs
STRINGdb), pseudotime (needs Monocle3), sctransform (marginal value).
"""
from pathlib import Path

from biopipen.core.proc import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import HdWGCNA as HdWGCNA_
from biopipen.core.testing import get_pipeline

DATA_DIR = Path(__file__).parent / "data"
DRIVE_ID = "1yxolklYrwFB9Snwr2Dp_W2eunBxaol4A"
# Name of the prepared object relative to the input object's directory,
# resolved by the process when `ref_srtobj` is not an existing file.
PREPARED_RDS = "Zhou_2020_control.prepared.rds"


class DownloadZhou2020(Proc):
    """Download the Zhou et al. 2020 control snRNA-seq dataset."""
    lang = config.lang.python
    input = "seed:var"
    input_data = [""]
    output = "outfile:file:Zhou_2020_control.rds"
    envs = {
        "drive_id": DRIVE_ID,
        "cache_path": str(DATA_DIR / "Zhou_2020_control.rds"),
    }
    script = """
        import shutil
        from pathlib import Path

        outfile = Path({{out.outfile | r}})
        cached = Path({{envs.cache_path | r}})
        if cached.exists():
            print(f"Using cached data at {cached} ...")
            shutil.copy2(cached, outfile)
        else:
            print("Downloading the Zhou et al. 2020 control dataset ...")
            import gdown
            try:
                gdown.download(id={{envs.drive_id | r}}, output=str(outfile))
            except Exception:
                gdown.download_folder(
                    id={{envs.drive_id | r}}, output=str(outfile.parent)
                )
            if not outfile.exists():
                raise RuntimeError(
                    "Failed to download the dataset. Please download "
                    "Zhou_2020_control.rds manually from "
                    "https://drive.google.com/file/d/1yxolklYrwFB9Snwr2Dp_W2eunBxaol4A/view "
                    f"and place it at {cached}"
                )
            cached.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(outfile, cached)
        print("Done")
    """


class PrepareSeurat(Proc):
    """Prepare the Seurat v5 object: UpdateSeuratObject + extra meta columns."""
    requires = DownloadZhou2020
    input = "srtfile:file"
    output = "outfile:file:{{in.srtfile | stem}}.prepared.rds"
    lang = config.lang.rscript
    script = """
        library(Seurat)
        srtfile = {{in.srtfile | r}}
        outfile = {{out.outfile | r}}
        srtobj = readRDS(srtfile)
        srtobj = UpdateSeuratObject(srtobj)
        # subset to a manageable size for local testing (keep all cell
        # types); the full ~20k-cell object is too heavy for WSL2
        set.seed(8525)
        cells = unlist(lapply(
            split(colnames(srtobj), srtobj$cell_type),
            function(cc) sample(cc, min(1000, length(cc)))
        ))
        srtobj = subset(srtobj, cells = cells)
        # the basic tutorial runs FindVariableFeatures + ScaleData during
        # preprocessing, which the uploaded object does not have
        if (length(VariableFeatures(srtobj)) == 0) {
            srtobj = FindVariableFeatures(srtobj, nfeatures = 2000)
        }
        if (!"ScaleData.RNA" %in% names(srtobj@commands)) {
            srtobj = ScaleData(srtobj)
        }
        meta = srtobj@meta.data
        # the pseudobulk consensus tutorial casts msex to a factor
        if ("msex" %in% colnames(meta) && !is.factor(meta$msex)) {
            meta$msex = as.factor(meta$msex)
        }
        # the network viz tutorial derives a `cluster` column from `annotation`
        if ("annotation" %in% colnames(meta) && !"cluster" %in% colnames(meta)) {
            meta$cluster = do.call(rbind, strsplit(as.character(meta$annotation), " "))[, 1]
        }
        # the module-trait tutorial uses `age_death` as a numeric trait, but
        # the uploaded object stores it as a character
        if ("age_death" %in% colnames(meta) && is.character(meta$age_death)) {
            meta$age_death = as.numeric(meta$age_death)
        }
        srtobj@meta.data = meta
        saveRDS(srtobj, outfile)
    """


# The INH metacell network envs from the basic tutorial, shared by the
# metacell-based test classes (basic, viz, enrichment, preservation, motifs,
# TF network).
BASIC_ENVS = {
    "SetupForWGCNA": {
        "wgcna_name": "tutorial",
        "gene_select": "fraction",
        "fraction": 0.05,
    },
    "MetacellsByGroups": {
        "group-by": ["cell_type", "Sample"],
        "ident-group": "cell_type",
        "k": 25,
        "max_shared": 10,
        "min_cells": 50,
        "reduction": "harmony",
    },
    "NormalizeMetacells": {},
    "ScaleMetacells": {"features": "r:VariableFeatures(srtobj)"},
    "RunPCAMetacells": {"features": "r:VariableFeatures(srtobj)"},
    "RunHarmonyMetacells": {"group-by-vars": "Sample"},
    "RunUMAPMetacells": {"reduction": "harmony", "dims": "1:15"},
    "SetDatExpr": {
        "group_name": "INH",
        "group-by": "cell_type",
        "assay": "RNA",
        "layer": "data",
    },
    "TestSoftPowers": {"networkType": "signed"},
    "ConstructNetwork": {"tom_name": "INH"},
    "ModuleEigengenes": {"group-by-vars": "Sample"},
    "ModuleConnectivity": {"group-by": "cell_type", "group_name": "INH"},
    "plots": {
        "Soft Powers": {"kind": "soft_powers"},
        "Dendrogram": {"kind": "dendrogram"},
    },
}


class TestHdWGCNAMetacells(HdWGCNA_):
    """Replicates the basic tutorial: harmony metacells + INH network.

    Tutorial: https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html
    """
    requires = PrepareSeurat
    envs = {
        **BASIC_ENVS,
        # the basic tutorial also resets the module names, computes the
        # hub-gene signature scores (ModuleExprScore, UCell), and plots the
        # module UMAP and the "Basic Visualization" section plots
        "ResetModuleNames": {"new_name": "INH-M"},
        "ModuleExprScore": {"n_genes": 25, "method": "UCell"},
        "plots": {
            **BASIC_ENVS["plots"],
            "Module UMAP": {
                "kind": "module_umap",
                "n_hubs": 10,
                "n_neighbors": 15,
                "min_dist": 0.1,
                "umap_plot_args": {
                    "edge-alpha": 0.25,
                    "sample_edges": True,
                    "edge_prop": 0.1,
                    "label_hubs": 2,
                    "keep_grey_edges": False,
                },
            },
            "hMEs": {
                "kind": "module_feature",
                "features": "hMEs",
                "order": True,
            },
            "Scores": {
                "kind": "module_feature",
                "features": "scores",
                "order": "shuffle",
                "ucell": True,
            },
            "Module Radar": {
                "kind": "module_radar",
                "group-by": "cluster",
                "barcodes": (
                    "r:rownames(srtobj@meta.data)[srtobj@meta.data$cell_type == 'INH']"
                ),
            },
            "Module Correlogram": {"kind": "module_correlogram"},
            "Module KMEs": {"kind": "kmes", "ncol": 5},
        },
    }


class TestHdWGCNAPseudobulk(HdWGCNA_):
    """Replicates the pseudobulk tutorial: VST pseudobulk network.

    Tutorial: https://smorabit.github.io/hdWGCNA/articles/pseudobulk.html
    """
    requires = PrepareSeurat
    envs = {
        "use_pseudobulk": True,
        "SetupForWGCNA": {
            "wgcna_name": "pseudobulk",
            "gene_select": "fraction",
            "fraction": 0.05,
        },
        "AggregatePseudobulk": {"replicate_col": "Sample", "group_col": "cell_type"},
        "NormalizeCounts": {},
        "SetDatExpr": {},
        "TestSoftPowers": {},
        "ConstructNetwork": {"tom_name": "pseudobulk", "mergeCutHeight": 0.15},
        "ModuleEigengenes": {},
        "ModuleConnectivity": {},
        "module_trait_corr": {
            "All": {"traits": ["msex", "braaksc"], "group-by": "cell_type"},
        },
        "plots": {
            "Soft Powers": {"kind": "soft_powers"},
            "Dendrogram": {"kind": "dendrogram"},
            "Module UMAP": {
                "kind": "module_umap",
                "n_hubs": 5,
                "n_neighbors": 10,
                "min_dist": 0.4,
                "spread": 3,
                "supervised": True,
                "target_weight": 0.3,
            },
        },
    }


class TestHdWGCNAConsensus(HdWGCNA_):
    """Replicates the consensus WGCNA tutorial.

    Tutorial: https://smorabit.github.io/hdWGCNA/articles/consensus_wgcna.html
    """
    requires = PrepareSeurat
    envs = {
        "use_pseudobulk": True,
        "use_consensus": True,
        "SetupForWGCNA": {
            "wgcna_name": "pseudobulk_consensus",
            "gene_select": "fraction",
            "fraction": 0.05,
        },
        "AggregatePseudobulk": {"replicate_col": "Sample", "group_col": "cell_type"},
        "NormalizeCounts": {},
        "SetMultiExpr": {"multi-group-by": "msex"},
        "TestSoftPowersConsensus": {},
        "ConstructNetwork": {"tom_name": "pseudobulk_consensus"},
        "ModuleEigengenes": {},
        "ModuleConnectivity": {},
        "plots": {
            "Soft Powers": {"kind": "soft_powers"},
            "Dendrogram": {"kind": "dendrogram"},
            "Module UMAP": {
                "kind": "module_umap",
                "n_hubs": 5,
                "n_neighbors": 5,
                "min_dist": 0.1,
                "spread": 2,
                "supervised": True,
                "target_weight": 0.05,
                "umap_plot_args": {
                    "edge-alpha": 0.5,
                    "sample_edges": True,
                    "keep_grey_edges": False,
                    "edge_prop": 0.075,
                    "label_hubs": 0,
                },
            },
            # the tutorial's comparison figure: the standard (non-consensus)
            # workflow re-run on the same data, both module-color tracks
            # plotted on the consensus dendrogram
            "Consensus vs Standard": {
                "kind": "consensus_compare",
                "standard_name": "pseudobulk_standard",
                "soft_power": 8,
                "main": "Sex consensus dendrogram",
            },
        },
    }


class TestHdWGCNADownstream(HdWGCNA_):
    """Replicates the DME and module-trait correlation tutorials.

    Both tutorials run on the metacell INH network from the basic tutorial.
    Tutorials:
    https://smorabit.github.io/hdWGCNA/articles/differential_MEs.html
    https://smorabit.github.io/hdWGCNA/articles/module_trait_correlation.html
    """
    requires = PrepareSeurat
    envs = {
        **BASIC_ENVS,
        "ModuleExprScore": {"n_genes": 25, "method": "UCell"},
        "dmes": {
            # one-vs-all DME analysis with FindAllDMEs (cell-type groups),
            # plotted as a volcano (the tutorial uses plot_labels and
            # show_cutoff to unclutter the plot)
            "FindAllDMEs": {
                "mode": "find_all",
                "group-by": "cell_type",
                "lollipop": False,
                "volcano": True,
                "plot_labels": False,
                "show_cutoff": False,
            },
            # module score DME analysis (FindAllDMEs on the hub-gene
            # signature scores from ModuleExprScore)
            "ModuleScores": {
                "mode": "find_all",
                "group-by": "cell_type",
                "features": "ModuleScores",
                "lollipop": True,
                "volcano": True,
            },
            # two-group DME analysis with FindDMEs (INH, female vs male)
            "INH by sex": {
                "mode": "find",
                "test-use": "wilcox",
                "barcodes1": (
                    "r:rownames(subset(srtobj@meta.data, cell_type == 'INH' & msex == 0))"
                ),
                "barcodes2": (
                    "r:rownames(subset(srtobj@meta.data, cell_type == 'INH' & msex != 0))"
                ),
                "lollipop": True,
                "volcano": True,
            },
        },
        "module_trait_corr": {
            "All": {
                "traits": [
                    "braaksc",
                    "pmi",
                    "msex",
                    "age_death",
                    "doublet_scores",
                    "nCount_RNA",
                    "nFeature_RNA",
                    "total_counts_mt",
                ],
                "group-by": "cell_type",
                "plot_args": {
                    "label": "fdr",
                    "label_symbol": "stars",
                    "text_size": 2,
                    "text_digits": 2,
                    "text_color": "white",
                    "high_color": "yellow",
                    "mid_color": "black",
                    "low_color": "purple",
                    "plot_max": 0.2,
                    "combine": True,
                },
            },
        },
    }


class TestHdWGCNAEnrichment(HdWGCNA_):
    """Replicates the enrichment tutorial (queries the Enrichr API).

    The tutorial uses the three Gene Ontology 2023 databases and the
    EnrichrBarPlot / EnrichrDotPlot visualizations (here replicated with
    enrichit RunEnrichment + VizEnrichment).
    Also runs a GSEA of the module kME ranks against the MSigDB C5:BP
    gene sets (fgsea + plotthis::VizGSEA).
    Tutorial: https://smorabit.github.io/hdWGCNA/articles/enrichment_analysis.html
    """
    requires = PrepareSeurat
    envs = {
        **BASIC_ENVS,
        "enrich": {
            "All": {
                "dbs": [
                    "GO_Biological_Process_2023",
                    "GO_Cellular_Component_2023",
                    "GO_Molecular_Function_2023",
                ],
                "enrich_plots": {
                    "Bar Plot": {
                        "db": "GO_Biological_Process_2023",
                        "plot_type": "bar",
                        "top_term": 10,
                        "ncol": 1,
                    },
                    "Dot Plot": {
                        "db": "GO_Cellular_Component_2023",
                        "plot_type": "dot",
                        "top_term": 10,
                        "ncol": 1,
                    },
                },
            },
        },
        "gsea": {
            "All": {
                "genesets": {
                    "species": "human",
                    "collection": "C5",
                    "subcollection": "BP",
                },
                "minSize": 10,
                "maxSize": 500,
                "top_term": 5,
                "gsea_plots": {
                    "Top Terms": {},
                },
            },
        },
    }


class TestHdWGCNAViz(HdWGCNA_):
    """Replicates the network visualizations tutorial.

    Tutorial: https://smorabit.github.io/hdWGCNA/articles/network_visualizations.html
    Note: the tutorial's custom ggraph/tidygraph network plots are not
    replicated (they are fully custom code, not hdWGCNA functions).
    """
    requires = PrepareSeurat
    envs = {
        **BASIC_ENVS,
        "plots": {
            **BASIC_ENVS["plots"],
            "Module Network": {"kind": "module_network"},
            "Hub Gene Network": {
                "kind": "hub_gene_network",
                "n_hubs": 3,
                "n_other": 5,
                "edge_prop": 0.75,
                "mods": "all",
            },
            "Module UMAP": {
                "kind": "module_umap",
                "n_hubs": 10,
                "n_neighbors": 15,
                "min_dist": 0.1,
                "umap_plot_args": {
                    "edge-alpha": 0.25,
                    "sample_edges": True,
                    "edge_prop": 0.1,
                    "label_hubs": 2,
                    "keep_grey_edges": False,
                },
            },
            "Module Feature": {
                "kind": "module_feature",
                "features": "hMEs",
                "order": True,
            },
            "Module Radar": {
                "kind": "module_radar",
                "group-by": "cluster",
                "barcodes": (
                    "r:rownames(srtobj@meta.data)[srtobj@meta.data$cell_type == 'INH']"
                ),
            },
            "Module Correlogram": {"kind": "module_correlogram"},
            "Module KMEs": {"kind": "kmes", "ncol": 5},
        },
    }


class TestHdWGCNAPreservation(HdWGCNA_):
    """Replicates the module preservation tutorial (query == reference).

    A cross-dataset preservation would need a second pre-built reference
    object; here we use the hdWGCNA-documented self-test pattern
    (`ref_srtobj: "self"`), preserving the modules against the query object
    itself. The tutorial's NetRep analysis (ModulePreservationNetRep, needs
    the NetRep package) and the module topology plots are also recovered.
    Tutorial: https://smorabit.github.io/hdWGCNA/articles/module_preservation.html
    """
    requires = PrepareSeurat
    envs = {
        **BASIC_ENVS,
        "ref_srtobj": "self",
        "module_preservations": {
            "Self": {
                "project_modules": True,
                "preserve": True,
                "plot": True,
                "statistics": "summary",
            },
            "NetRep": {
                "project_modules": False,
                "preserve": False,
                "plot": False,
                "netrep": {
                    "args": {"n_permutations": 100, "n_threads": "r:ncores"},
                    "topology_heatmap": {
                        "Heatmap": {
                            "mod": (
                                "r:setdiff(unique(GetModules(srtobj)$module), 'grey')[1]"
                            ),
                            "matrix": "Cor",
                            "order_by": "degree",
                            "high_color": "red",
                            "plot_max": 0.75,
                        },
                    },
                    "topology_barplot": {
                        "Barplot": {
                            "mod": (
                                "r:setdiff(unique(GetModules(srtobj)$module), 'grey')[1]"
                            ),
                            "features": "weighted_degree",
                        },
                    },
                },
            },
        },
    }


class TestHdWGCNAMotifs(HdWGCNA_):
    """Replicates the motif analysis tutorial (needs Bioc genome packages).

    Tutorial: https://smorabit.github.io/hdWGCNA/articles/motif_analysis.html
    Note: the tutorial's `MotifTargetScore` and `ModuleTFNetwork` sections
    are not replicated because those functions are not exported by hdWGCNA
    0.4.12 (the tutorial is marked Deprecated). `MotifOverlapBarPlot` is also
    broken in 0.4.12 and is replicated in-house in HdWGCNA-motifs.R.
    """
    requires = PrepareSeurat
    envs = {
        **BASIC_ENVS,
        "tf_network": {
            "motif_scan": {
                "species_genome": "hg38",
                "ensdb_package": "EnsDb.Hsapiens.v86",
            },
        },
        "motifs": {
            "overlap_bar_plot": {},
        },
    }


class TestHdWGCNA_TFNetwork(HdWGCNA_):
    """Replicates the TF network tutorial (needs xgboost, Bioc packages).

    Tutorial: https://smorabit.github.io/hdWGCNA/articles/tf_network.html
    """
    requires = PrepareSeurat
    envs = {
        **BASIC_ENVS,
        "tf_network": {
            "motif_scan": {
                "species_genome": "hg38",
                "ensdb_package": "EnsDb.Hsapiens.v86",
            },
            "construct": {
                "model_params": {
                    "objective": "reg:squarederror",
                    "max_depth": 1,
                    "eta": 0.1,
                    "alpha": 0.5,
                },
            },
            "assign": {
                "strategy": "A",
                "reg_thresh": 0.01,
                "n_tfs": 10,
            },
            # without `target_type`, the scores are computed for both
            # "positive" and "negative" (required by FindDifferentialRegulons)
            "regulon_scores": {"ncores": 8},
            "plots": {
                "Network": {
                    "selected_tfs": "r:c('RUNX2', 'RXRA', 'TCF4')",
                    "label_TFs": 0,
                    "depth": 1,
                },
            },
            "regulon_bar_plots": {
                "Bar": {"selected_tf": "RUNX2"},
            },
            "differential_regulons": {
                # female vs male in the inhibitory neurons, as in the
                # tutorial
                "INH by sex": {
                    "barcodes1": (
                        "r:rownames(subset(srtobj@meta.data, cell_type == 'INH' & msex == 0))"
                    ),
                    "barcodes2": (
                        "r:rownames(subset(srtobj@meta.data, cell_type == 'INH' & msex != 0))"
                    ),
                },
            },
            "enrich_regulons": {"dbs": ["GO_Biological_Process_2021"]},
            "regulatory_heatmap": {"feature": "delta", "dendrogram": False},
            "regulatory_network": {"feature": "delta", "cutoff": 0.5, "max_val": 1.5},
        },
        # the module UMAP is required by ModuleRegulatoryNetworkPlot
        "plots": {
            **BASIC_ENVS["plots"],
            "Module UMAP": {"kind": "module_umap"},
        },
    }


def pipeline():
    return get_pipeline(__file__, enable_report=True).set_starts(DownloadZhou2020)


def testing(pipen):
    import json
    import struct

    outdir = Path(pipen.outdir)

    def jobdir(cls):
        return outdir / cls.__name__ / "Zhou_2020_control.prepared.hdwgcn"

    def assert_png_size(path, min_width=800, min_height=600):
        """Multi-panel plots must not be squeezed into degenerate devices."""
        with open(path, "rb") as fh:
            assert fh.read(8) == b"\x89PNG\r\n\x1a\n", f"Not a PNG: {path}"
            fh.seek(16)
            width, height = struct.unpack(">II", fh.read(8))
        assert width >= min_width and height >= min_height, (
            f"PNG too small ({width}x{height}): {path}"
        )

    def assert_no_abs_paths(cls):
        """Report descriptions must not leak the server (output) paths.

        The `src` fields are absolute by convention (so do all processes)
        and are relativized when the report is built, but the descriptions
        are rendered verbatim.
        """
        report = json.loads((outdir / cls.__name__ / "report.json").read_text())

        def walk(node):
            if isinstance(node, dict):
                if (
                    node.get("kind") == "descr"
                    and outdir.as_posix() in str(node.get("content", ""))
                ):
                    raise AssertionError(
                        f"Absolute path in report descr: {str(node.get('content'))[:120]}"
                    )
                for v in node.values():
                    walk(v)
            elif isinstance(node, list):
                for v in node:
                    walk(v)

        walk(report)

    # common outputs of every network construction
    for cls in [
        TestHdWGCNAMetacells,
        TestHdWGCNAPseudobulk,
        TestHdWGCNAConsensus,
        TestHdWGCNADownstream,
        TestHdWGCNAEnrichment,
        TestHdWGCNAViz,
        TestHdWGCNAPreservation,
        TestHdWGCNAMotifs,
        TestHdWGCNA_TFNetwork,
    ]:
        jdir = jobdir(cls)
        assert (jdir / "Zhou_2020_control.prepared.hdwgcn.qs").exists()
        mod_table = (jdir / "tables" / "modules.tsv").read_text()
        assert "gene_name" in mod_table
        assert "module" in mod_table
        assert "color" in mod_table
        assert "kME_" in mod_table
        assert (jdir / "tables" / "hub_genes.tsv").exists()
        assert (jdir / "tables" / "power_table.tsv").exists()
        assert (jdir / "tables" / "MEs.tsv").exists()
        assert_no_abs_paths(cls)

    # basic tutorial: module names reset, module UMAP and the "Basic
    # Visualization" plots (module features, radar, correlogram, KMEs)
    jdir = jobdir(TestHdWGCNAMetacells)
    assert "INH-M" in (jdir / "tables" / "modules.tsv").read_text()
    assert (jdir / "plots" / "Soft-Powers" / "Soft-Powers.sft.png").exists()
    assert (jdir / "plots" / "Soft-Powers" / "Soft-Powers.fit.png").exists()
    assert (jdir / "plots" / "Dendrogram" / "Dendrogram.png").exists()
    assert (jdir / "plots" / "Module-UMAP" / "Module-UMAP.png").exists()
    assert (jdir / "plots" / "Module-Features" / "hMEs.png").exists()
    assert (jdir / "plots" / "Module-Features" / "Scores.png").exists()
    assert (jdir / "plots" / "Module-Radars" / "Module-Radar.png").exists()
    assert (jdir / "plots" / "Module-Correlograms" / "Module-Correlogram.png").exists()
    assert (jdir / "plots" / "Module-KMEs" / "Module-KMEs.png").exists()

    # pseudobulk tutorial
    jdir = jobdir(TestHdWGCNAPseudobulk)
    assert (jdir / "plots" / "Soft-Powers" / "Soft-Powers.sft.png").exists()
    assert (jdir / "plots" / "Dendrogram" / "Dendrogram.png").exists()
    assert (jdir / "plots" / "Module-UMAP" / "Module-UMAP.png").exists()
    mtc_dir = jdir / "tables" / "module_trait_corr" / "Module-Trait-Correlation"
    assert (mtc_dir / "All.cor.tsv").exists()
    assert (mtc_dir / "All.png").exists()

    # consensus tutorial
    jdir = jobdir(TestHdWGCNAConsensus)
    assert (jdir / "plots" / "Soft-Powers" / "Soft-Powers.sft.png").exists()
    assert (jdir / "plots" / "Dendrogram" / "Dendrogram.png").exists()
    assert (jdir / "plots" / "Module-UMAP" / "Module-UMAP.png").exists()
    assert (jdir / "tables" / "power_table.tsv").exists()
    # the consensus tutorial's per-set soft powers (each plot type shows both
    # sets side-by-side) and the consensus-vs-standard comparison dendrogram
    sp_dir = jdir / "plots" / "Soft-Powers"
    for suffix in ("sft", "fit", "median.k", "max.k"):
        assert_png_size(sp_dir / f"Soft-Powers.{suffix}.png")
    assert (jdir / "plots" / "Dendrogram" / "Consensus-vs-Standard.png").exists()
    # the standard workflow dendrogram (tutorial figure 3)
    assert (jdir / "plots" / "Dendrogram" / "Consensus-vs-Standard.standard.png").exists()

    # DME + module-trait correlation
    jdir = jobdir(TestHdWGCNADownstream)
    dmes_dir = jdir / "tables" / "dmes" / "DMEs"
    assert (dmes_dir / "FindAllDMEs.tsv").exists()
    assert (dmes_dir / "FindAllDMEs.volcano.png").exists()
    assert (dmes_dir / "ModuleScores.tsv").exists()
    assert_png_size(dmes_dir / "ModuleScores.lollipop.png")
    assert (dmes_dir / "ModuleScores.volcano.png").exists()
    assert (dmes_dir / "INH-by-sex.tsv").exists()
    assert (dmes_dir / "INH-by-sex.lollipop.png").exists()
    assert (dmes_dir / "INH-by-sex.volcano.png").exists()
    mtc_dir = jdir / "tables" / "module_trait_corr" / "Module-Trait-Correlation"
    assert (mtc_dir / "All.cor.tsv").exists()
    assert (mtc_dir / "All.pval.tsv").exists()
    assert (mtc_dir / "All.fdr.tsv").exists()
    assert (mtc_dir / "All.png").exists()

    # enrichment
    jdir = jobdir(TestHdWGCNAEnrichment)
    enrich_dir = jdir / "tables" / "enrich" / "Enrichment"
    tsvs = list(enrich_dir.glob("All.*.tsv"))
    assert tsvs, f"No enrichment tables in {enrich_dir}"
    assert list((jdir / "plots" / "Enrich").glob("All.*.png")), "No enrichment plots"

    # GSEA of the module kME ranks
    gsea_dir = jdir / "tables" / "gsea" / "GSEA"
    gsea_tsvs = list(gsea_dir.glob("All.*.tsv"))
    assert gsea_tsvs, f"No GSEA tables in {gsea_dir}"
    gsea_plots = jdir / "plots" / "GSEA"
    assert list(gsea_plots.glob("All.*.summary.png")), "No GSEA summary plots"
    assert list(gsea_plots.glob("All.*.Top-Terms.png")), "No GSEA running-score plots"
    # the summary dot plot with per-term line plots keeps its natural size
    # (from the plotthis attrs, ~6.5in wide); it must not have fallen back
    # to the default 800x600 device box
    with open(gsea_plots / gsea_plots.glob("All.*.summary.png").__next__().name, "rb") as fh:
        assert fh.read(8) == b"\x89PNG\r\n\x1a\n", "Not a PNG"
        fh.seek(16)
        gsea_width, gsea_height = struct.unpack(">II", fh.read(8))
    assert gsea_width > 800 or gsea_height > 600, (
        f"GSEA summary PNG is the default device box ({gsea_width}x{gsea_height})"
    )

    # network visualizations
    jdir = jobdir(TestHdWGCNAViz)
    # ModuleNetworkPlot writes one PDF per module; the report embeds a
    # single merged PDF instead
    assert (jdir / "plots" / "ModuleNetworks.Module-Network" / "module_networks.pdf").exists()
    report = json.loads((outdir / TestHdWGCNAViz.__name__ / "report.json").read_text())
    net_items = [
        item
        for sec in report.get("Module Networks", {}).values()
        for tab in sec.get("#", {}).get("tabs", [])
        for item in tab.get("contents", [])
    ]
    assert any(
        item.get("kind") == "pdf" and item.get("src", "").endswith("module_networks.pdf")
        for item in net_items
    ), f"No merged module-network PDF in the report: {[i.get('kind') for i in net_items]}"
    assert (jdir / "plots" / "Hub-Gene-Networks" / "Hub-Gene-Network.png").exists()
    assert (jdir / "plots" / "Module-Features" / "Module-Feature.png").exists()
    assert (jdir / "plots" / "Module-Radars" / "Module-Radar.png").exists()
    assert (
        jdir / "plots" / "Module-Correlograms" / "Module-Correlogram.png"
    ).exists()
    assert (jdir / "plots" / "Module-UMAP" / "Module-UMAP.png").exists()
    # multi-panel plots must get sized devices, not the fixed 800x600 box
    assert_png_size(jdir / "plots" / "Module-KMEs" / "Module-KMEs.png")
    assert_png_size(jdir / "plots" / "Module-Radars" / "Module-Radar.png")
    assert_png_size(jdir / "plots" / "Module-Features" / "Module-Feature.png")

    # module preservation
    jdir = jobdir(TestHdWGCNAPreservation)
    pres_dir = jdir / "tables" / "preservation" / "Module-Preservation"
    assert (pres_dir / "Self.Z.tsv").exists()
    assert (pres_dir / "Self.obs.tsv").exists()
    assert_png_size(pres_dir / "Self.png")
    # NetRep topology plots
    topo_dir = jdir / "plots" / "Topology" / "Module-Preservation"
    assert (topo_dir / "NetRep-Heatmap.png").exists()
    assert (topo_dir / "NetRep-Barplot.png").exists()
    # every preservation tab must have real content: a pure-netrep case
    # (no preserve/plot of its own) must not add an empty tab
    report = json.loads((outdir / TestHdWGCNAPreservation.__name__ / "report.json").read_text())
    for h2, sec in report.get("Module Preservation", {}).items():
        for tab in sec.get("#", {}).get("tabs", []):
            # tabs must be dicts; a list (from passing a list of tabs as a
            # single add2 argument) breaks the report render
            assert isinstance(tab, dict) and "contents" in tab, (
                f"Malformed preservation tab: {h2}/{tab}"
            )
            items = tab.get("contents", [])
            assert any(item.get("kind") != "descr" for item in items), (
                f"Empty preservation tab: {h2}/{tab.get('name')}"
            )

    # motifs
    jdir = jobdir(TestHdWGCNAMotifs)
    assert (jdir / "tables" / "motif_overlap.tsv").exists()
    # the per-module PDFs are merged into one; the report embeds that single
    # file instead of one PDF per module
    assert (jdir / "plots" / "MotifOverlap" / "motif_overlaps.pdf").exists()
    report = json.loads((outdir / TestHdWGCNAMotifs.__name__ / "report.json").read_text())
    pdf_items = [
        item
        for sec in report.get("Motifs", {}).values()
        for tab in sec.get("#", {}).get("tabs", [])
        for item in tab.get("contents", [])
        if item.get("kind") == "pdf"
    ]
    assert len(pdf_items) == 1, f"Expected one merged motif PDF, got: {pdf_items}"
    assert pdf_items[0]["src"].endswith("motif_overlaps.pdf"), pdf_items[0]

    # TF network
    jdir = jobdir(TestHdWGCNA_TFNetwork)
    assert (jdir / "plots" / "TF" / "TF-Network" / "Network.png").exists()
    assert (jdir / "plots" / "TF" / "TF-Network" / "Bar.png").exists()
    tf_tables = jdir / "tables" / "tf_network"
    assert (tf_tables / "TF-Network" / "INH-by-sex.tsv").exists()
    assert (tf_tables / "regulon_enrichment.tsv").exists()
    assert (
        jdir / "plots" / "TF" / "Module-Regulatory-Heatmap" / "Regulatory-Heatmap.png"
    ).exists()
    assert (
        jdir / "plots" / "TF" / "Module-Regulatory-Network" / "Regulatory-Network.png"
    ).exists()


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
