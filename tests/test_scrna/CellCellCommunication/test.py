from biopipen.core.proc import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import (
    CellCellCommunication as CellCellCommunication_,
    CellCellCommunicationPlots as CellCellCommunicationPlots_,
    AnnData2Seurat as AnnData2Seurat_,
)
from biopipen.core.testing import get_pipeline


class PrepareAnnData(Proc):
    lang = config.lang.python
    input = "var"
    input_data = ["toy"]
    output = "outfile:file:{{in.var}}.h5ad"
    script = """
        import scanpy as sc

        adata = sc.datasets.pbmc68k_reduced()
        adata.write_h5ad({{out.outfile | repr}})
    """


class CellCellCommunicationAnndata(CellCellCommunication_):
    requires = PrepareAnnData
    envs = {"groupby": "bulk_labels"}


class CellCellCommunicationAnndataPlots(CellCellCommunicationPlots_):
    requires = CellCellCommunicationAnndata
    envs = {
        "cases": {
            "Heatmap": {
                "kind": "heatmap",
                "section": "Heatmap",
            },
            "HeatmapB": {
                "kind": "heatmap",
                "option": "B",
                "n_top_ints": 10,
                "section": "Heatmap",
            },
            "HeatmapCellPhoneDB": {
                "kind": "heatmap",
                "option": "CellPhoneDB",
                "section": "Heatmap",
            },
            "DotPlot": {
                "kind": "dotplot",
                "section": "DotPlot",
            },
            "DotPlotB": {
                "kind": "dotplot",
                "option": "B",
                "section": "DotPlot",
            },
            "DotPlotLiana": {
                "kind": "dotplot",
                "option": "Liana",
                "section": "DotPlot",
            },
            "Network": {
                "kind": "network",
                "section": "Network",
            },
            "NetworkB": {
                "kind": "network",
                "option": "B",
                "section": "Network",
            },
            "Circos": {
                "kind": "circos",
                "section": "Circos",
            },
            "CircosB": {
                "kind": "circos",
                "option": "B",
                "section": "Circos",
            },
            "CircosC": {
                "kind": "circos",
                "option": "C",
                "section": "Circos",
            },
            "Arrow": {
                "kind": "arrow",
                "section": "Arrow",
                "cell_types": ["CD34+", "CD19+ B"],
            },
            "ArrowB": {
                "kind": "arrow",
                "section": "Arrow",
                "option": "B",
                "cell_types": ["CD34+", "CD19+ B"],
            },
            "Sigmoid": {
                "kind": "sigmoid",
                "section": "Sigmoid",
            },
        }
    }


class AnnData2Seurat(AnnData2Seurat_):
    requires = PrepareAnnData


class CellCellCommunicationSeurat(CellCellCommunication_):
    requires = AnnData2Seurat
    envs = {"groupby": "bulk_labels"}


class CellCellCommunicationSeuratPlots(CellCellCommunicationPlots_):
    requires = CellCellCommunicationSeurat
    envs = {
        "cases": {
            "Heatmap": {
                "kind": "heatmap",
                "section": "Heatmap",
            },
            "HeatmapB": {
                "kind": "heatmap",
                "option": "B",
                "n_top_ints": 10,
                "section": "Heatmap",
            },
            "HeatmapCellPhoneDB": {
                "kind": "heatmap",
                "option": "CellPhoneDB",
                "section": "Heatmap",
            },
            "DotPlot": {
                "kind": "dotplot",
                "section": "DotPlot",
            },
            "DotPlotB": {
                "kind": "dotplot",
                "option": "B",
                "section": "DotPlot",
            },
            "DotPlotLiana": {
                "kind": "dotplot",
                "option": "Liana",
                "section": "DotPlot",
            },
            "Network": {
                "kind": "network",
                "section": "Network",
            },
            "NetworkB": {
                "kind": "network",
                "option": "B",
                "section": "Network",
            },
            "Circos": {
                "kind": "circos",
                "section": "Circos",
            },
            "CircosB": {
                "kind": "circos",
                "option": "B",
                "section": "Circos",
            },
            "CircosC": {
                "kind": "circos",
                "option": "C",
                "section": "Circos",
            },
            "Arrow": {
                "kind": "arrow",
                "section": "Arrow",
                "cell_types": ["CD34+", "CD19+ B"],
            },
            "ArrowB": {
                "kind": "arrow",
                "section": "Arrow",
                "option": "B",
                "cell_types": ["CD34+", "CD19+ B"],
            },
            "Sigmoid": {
                "kind": "sigmoid",
                "section": "Sigmoid",
            },
        }
    }


def pipeline():
    return get_pipeline(__file__).set_starts(PrepareAnnData)


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
