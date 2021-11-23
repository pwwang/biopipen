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


# class FGSEA(Proc):


class Enrichr(Proc):
    """Gene set enrichment analysis using Enrichr

    Need `devtools::install_github("wjawaid/enrichR")`

    Input:
        infile: The gene list file.
            You can specify whether this file has header and the index (0-based)
            of the columns where the genes are present

    Output:
        outdir: The output directory

    Envs:
        inheader: Whether the input file has header (first row skipped)
        incol: Which column has the genes (0-based index)
        dbs: The databases to enrich against.
            See https://maayanlab.cloud/Enrichr/#libraries for all available
            databases/libaries
    """
    input = "infile:file"
    output = "outdir:dir:{{in.infile | stem}}.enrichr"
    lang = config.lang.rscript
    envs = {
        "inheader": False,
        "incol": 0,
        "genename": "symbol",
        "dbs": ["KEGG_2021_Human"],
    }
    script = "file://../scripts/gsea/Enrichr.R"
    plugin_opts = {
        "report": "file://../reports/gsea/Enrichr.svelte"
    }
