"""Tools to analyze single-cell RNA"""

from ..core.proc import Proc
from ..core.config import config

class SeuratLoading(Proc):
    """Seurat - Loading data

    Input:
        metafile: The metadata of the samples
            A tab-delimited file
            Two columns are required:
            - `Sample` to specify the sample names.
            - `RNADir` to assign the path of the data to the samples
              The path will be read by `Read10X()` from `Seurat`

    Output:
        rdsfile: The RDS file with a list of Seurat object
    """
    input = "metafile:file"
    output = "rdsfile:file:{{in.metafile | stem}}.seurat.RDS"
    lang = config.lang.rscript
    script = "file://../scripts/scrna/SeuratLoading.R"

class GeneExpressionInvistigation(Proc):
    """Investigation of expressions of genes of interest"""
    input = "srtobj:file, genefile:file, groupfile:file"
    output = "outdir:dir:{{in.groupfile | stem0}}.exprs"
    lang = config.lang.rscript
    envs = { "cases": {} }
    order = 4
    script = "file://../scripts/scrna/GeneExpressionInvistigation.R"
    plugin_opts = {
        "report": "file://../reports/scrna/GeneExpressionInvistigation.svelte"
    }
