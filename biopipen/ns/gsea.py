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
            `doc.string`: Documentation string used as a prefix
            to name result files. If not provided, will use
            `envs['doc.string']`

    Output:
        outdir: The output directory

    Envs:
        inopts: The options for `read.table()` to read the input file
            If `rds` will use `readRDS()`
        metaopts: The options for `read.table()` to read the meta file
        clscol: The column of the metafile determining the classes
        doc.string: Documentation string used as a prefix to name result files
            Other configs passed to `GSEA()` directly

    Requires:
        - name: GSEA-MSigDB/GSEA_R
          check: |
            {{proc.lang}} <(echo "library(GSEA)")
    """
    input = "infile:file, metafile:file, gmtfile:file, configfile:file"
    output = "outdir:dir:{{in.infile | stem}}.gsea"
    lang = config.lang.rscript
    envs = {
        "inopts": {"header": True, "row.names": -1},
        "metaopts": {"header": True, "row.names": -1},
        "clscol": None,
        "doc.string": "gsea_result",
    }
    script = "file://../scripts/gsea/GSEA.R"
    plugin_opts = {
        "report": "file://../reports/gsea/GSEA.svelte"
    }


class PreRank(Proc):
    """PreRank the genes for GSEA analysis

    Input:
        infile: The expression file.
            Either a tab-delimited matrix or an RDS file (on envs.inopts)
        metafile: The meta data file, determining the class of the samples
            Two columns are required
            Sample: The unique sample id for each sample
            `[Group]`: The groups/classes of the samples
        configfile: The configuration file in TOML format to specify some envs.
            `clscol`: If not provided, will use `envs.clscol`
            `classes`: Defines pos and neg labels. If not provided, use will
            `envs.classes`.

    Output:
        outfile: The rank file with 1st column the genes, and the rest the
            ranks for different class pairs provided by `envs.classes` or
            `in.configfile`

    Envs:
        inopts: Options for `read.table()` to read `in.infile`
        metaopts: Options for `read.table()` to read `in.metafile`
        method: The method to do the preranking.
            Supported: `s2n(signal_to_noise)`, `abs_s2n(abs_signal_to_noise)`,
            `t_test`, `ratio_of_classes`, `diff_of_classes` and
            `log2_ratio_of_classes`.
        clscol: The column of metafile specifying the classes of the samples
        classes: The classes to specify the pos and neg labels.
            It could be a pair of labels (e.g. `["CASE", "CNTRL"]`), where
            the first one is pos and second is neg. Or you can have multiple
            pairs of labels (e.g. `[["CASE1", "CNTRL"], ["CASE2", "CNTRL"]]`)
    """
    input = "infile:file, metafile:file, configfile:file"
    output = "outfile:file:{{in.infile | stem}}.rank"
    lang = config.lang.rscript
    envs = {
        "inopts": {"header": True, "row.names": -1},
        "metaopts": {"header": True, "row.names": -1},
        "method": "s2n",
        "clscol": None,
        "classes": None,
    }
    script = "file://../scripts/gsea/PreRank.R"


class FGSEA(Proc):
    """Gene set enrichment analysis using `fgsea`

    Need `devtools::install_github("ctlab/fgsea")`

    Input:
        infile: The expression file.
            Either a tab-delimited matrix or an RDS file (on envs.inopts)
        metafile: The meta data file, determining the class of the samples
            Two columns are required
            Sample: The unique sample id for each sample
            `[Group]`: The groups/classes of the samples
        gmtfile: The GMT file of reference gene sets
        configfile: The configuration file in TOML format to specify some envs.
            `clscol`: If not provided, will use `envs.clscol`
            `classes`: Defines pos and neg labels. If not provided, use will
            `envs.classes`.

    Output:
        outdir: The output directory

    Envs:
        inopts: The options for `read.table()` to read the input file
            If `rds` will use `readRDS()`
        metaopts: The options for `read.table()` to read the meta file
        method: The method to do the preranking.
            Supported: `s2n(signal_to_noise)`, `abs_s2n(abs_signal_to_noise)`,
            `t_test`, `ratio_of_classes`, `diff_of_classes` and
            `log2_ratio_of_classes`.
        clscol: The column of metafile specifying the classes of the samples
        classes: The classes to specify the pos and neg labels.
            It could be a pair of labels (e.g. `["CASE", "CNTRL"]`), where
            the first one is pos and second is neg. Or you can have multiple
            pairs of labels (e.g. `[["CASE1", "CNTRL"], ["CASE2", "CNTRL"]]`)
        top: Do gsea table and enrich plot for top N pathways. If it is < 1,
            will apply it to `padj`
        `<rest>`: Rest arguments for `fgsea()`

    Requires:
        - name: bioconductor-fgsea
          check: |
            {{proc.lang}} -e "library(fgsea)"
    """
    input = "infile:file, metafile:file, gmtfile:file, configfile:file"
    output = "outdir:dir:{{in.infile | stem}}.fgsea"
    lang = config.lang.rscript
    envs = {
        "inopts": {"header": True, "row.names": -1},
        "metaopts": {"header": True, "row.names": -1},
        "method": "s2n",
        "clscol": None,
        "classes": None,
        "top": 20,
        "ncores": config.misc.ncores,
        "minSize": 10,
        "maxSize": 100,
        "eps": 0,
    }
    script = "file://../scripts/gsea/FGSEA.R"
    plugin_opts = {
        "report": "file://../reports/gsea/FGSEA.svelte"
    }

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
        inopts: Options for `read.table()` to read `in.infile`
        genecol: Which column has the genes (0-based index or column name)
        dbs: The databases to enrich against.
            See https://maayanlab.cloud/Enrichr/#libraries for all available
            databases/libaries
    """
    input = "infile:file"
    output = "outdir:dir:{{in.infile | stem}}.enrichr"
    lang = config.lang.rscript
    envs = {
        "inopts": {},
        "genecol": 0,
        "genename": "symbol",
        "dbs": ["KEGG_2021_Human"],
    }
    script = "file://../scripts/gsea/Enrichr.R"
    plugin_opts = {
        "report": "file://../reports/gsea/Enrichr.svelte"
    }
