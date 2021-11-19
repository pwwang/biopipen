"""Gene set enrichment analysis"""

from ..core.proc import Proc
from ..core.config import config

class GSEA(Proc):
    """Gene set enrichment analysis

    Need `devtools::install_github("GSEA-MSigDB/GSEA_R")`

    Input:
        infile: The expression file.
            Either a tab-delimited matrix or an RDS file (on envs.infmt)
        metafile: The meta data file, determining the class of the samples
            Two columns are required
            Sample: The unique sample id for each sample
            `[Group]`: The groups/classes of the samples
        gmtfile: The GMT file of reference gene sets
        configfile: The configuration file in TOML format to specify some envs.
            `clscol`: If not provided, will use `envs.clscol`
            `doc.string` or `doc_string`: Documentation string used as a prefix
            to name result files. If not provided, will use
            `envs.doc_string` or `envs['doc_string']`

    Output:
        outdir: The output directory

    Envs:
        infmt: The format of the input file
            Either a tab-delimited matrxi file or an RDS file
        clscol: The column of the metafile determining the classes
        doc_string: Documentation string used as a prefix to name result files
            Other configs passed to `GSEA()` directly
    """
    input = "infile:file, metafile:file, gmtfile:file, configfile:file"
    output = "outdir:dir:{{in.infile | stem}}.gsea"
    lang = config.lang.rscript
    envs = {
        "infmt": "matrix",
        "clscol": None,
        "doc_string": "gsea_result",
    }
    script = "file://../scripts/gsea/GSEA.R"
    plugin_opts = {
        "report": "file://../reports/gsea/GSEA.svelte"
    }
