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
        adata.write_h5ad({{out.outfile | quote}})
    """


class CellCellCommunicationAnndata(CellCellCommunication_):
    requires = PrepareAnnData
    envs = {"groupby": "bulk_labels"}


class CellCellCommunicationAnndataPlots(CellCellCommunicationPlots_):
    requires = CellCellCommunicationAnndata
    envs = {
        "cases": {
            "Heatmap::Heatmap": {"plot_type": "heatmap"},
            "Heatmap::HeatmapB": {"plot_type": "heatmap", "method": "interaction"},
            "DotPlot::DotPlot": {"plot_type": "dot"},
            "DotPlot::DotPlotB": {"plot_type": "dot", "method": "interaction"},
            "Network::Network": {"plot_type": "network"},
            "Circos::Circos": {"plot_type": "circos"},
        }
    }


class AnnData2Seurat(AnnData2Seurat_):
    requires = PrepareAnnData
    envs = {"ident": "bulk_labels"}


class CellCellCommunicationSeurat(CellCellCommunication_):
    requires = AnnData2Seurat
    # envs = {"groupby": "bulk_labels"}


class CellCellCommunicationSeuratPlots(CellCellCommunicationPlots_):
    requires = CellCellCommunicationSeurat
    envs = {
        "cases": {
            "Heatmap::Heatmap": {"plot_type": "heatmap"},
            "Heatmap::HeatmapB": {"plot_type": "heatmap", "method": "interaction"},
            "DotPlot::DotPlot": {"plot_type": "dot"},
            "DotPlot::DotPlotB": {"plot_type": "dot", "method": "interaction"},
            "Network::Network": {"plot_type": "network"},
            "Circos::Circos": {"plot_type": "circos"},
        }
    }


def pipeline():
    return (
        get_pipeline(__file__)
        # get_pipeline(__file__, enable_report=True)
        .set_starts(PrepareAnnData)
    )


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
