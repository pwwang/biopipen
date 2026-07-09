from biopipen.core.proc import Proc
from biopipen.core.config import config
import pandas as pd
from pathlib import Path
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
    envs = {
        "groupby": "bulk_labels",
        "cases": {"DEFAULT": {}, "bulk_labels": {"groupby": "bulk_labels"}},
    }


class CellCellCommunicationAnndataSplitBy(CellCellCommunication_):
    requires = PrepareAnnData
    envs = {"groupby": "bulk_labels", "split_by": "phase"}


class CellCellCommunicationAnndataSubsetPython(CellCellCommunication_):
    requires = PrepareAnnData
    envs = {
        "groupby": "bulk_labels",
        "subset_using": "python",
        "subset": "louvain == '1'",
    }


class CellCellCommunicationAnndataPlots(CellCellCommunicationPlots_):
    requires = CellCellCommunicationAnndata
    envs = {
        "ligand_expr": "ligand_trimean",
        "receptor_expr": "receptor_trimean",
        "cases": {
            "Heatmap::Heatmap": {"plot_type": "heatmap", "subset": "Case == 'DEFAULT'"},
            "Heatmap::HeatmapB": {"plot_type": "heatmap", "method": "interaction"},
            "Heatmap::LinkedHeatmap": {"plot_type": "linkedheatmap"},
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
        "ligand_expr": "ligand_trimean",
        "receptor_expr": "receptor_trimean",
        "cases": {
            "Heatmap::Heatmap": {"plot_type": "heatmap"},
            "Heatmap::HeatmapB": {"plot_type": "heatmap", "method": "interaction"},
            "Heatmap::LinkedHeatmap": {"plot_type": "linkedheatmap"},
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
    outdir = Path(pipen.outdir)
    split_res = pd.read_csv(
        outdir / "CellCellCommunicationAnndataSplitBy" / "toy-ccc.txt",
        sep="\t",
    )
    assert "phase" in split_res.columns
    assert split_res["phase"].nunique() >= 2

    subset_res = pd.read_csv(
        outdir / "CellCellCommunicationAnndataSubsetPython" / "toy-ccc.txt",
        sep="\t",
    )
    assert not subset_res.empty


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
