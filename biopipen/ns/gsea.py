"""Gene set enrichment analysis"""
from pipen.utils import mark
from ..core.proc import Proc
from ..core.config import config


@mark(deprecated='[{proc.name}] is deprecated, use `FGSEA` instead.')
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
        doc_string: Documentation string used as a prefix to name result files
            Other configs passed to `GSEA()` directly

    Requires:
        GSEA-MSigDB/GSEA_R:
            - check: {{proc.lang}} <(echo "library(GSEA)")
    """

    input = "infile:file, metafile:file, gmtfile:file, configfile:file"
    output = "outdir:dir:{{in.infile | stem}}.gsea"
    lang = config.lang.rscript
    envs = {
        "inopts": {"header": True, "row.names": -1},
        "metaopts": {"header": True, "row.names": -1},
        "clscol": None,
        "doc_string": "gsea_result",
    }
    script = "file://../scripts/gsea/GSEA.R"
    plugin_opts = {"report": "file://../reports/gsea/GSEA.svelte"}


@mark(deprecated='[{proc.name}] is deprecated, use `FGSEA` directly.')
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

    Input:
        infile: The expression file (genes x samples).
            Either a tab-delimited file.
        metafile: The meta data file, determining the class of the samples
            Two columns are required. If column `Sample` is found, it will be used
            as the samples; otherwise the first column should be the samples.
            The other column should be the group/class of the samples, whose
            name is specified by `envs.clscol`.

    Output:
        outdir: The output directory containing the results, including
            the table and plots.

    Envs:
        ncores (type=int): Number of cores for parallelization
            Passed to `nproc` of `fgseaMultilevel()`.
        case: The case label for the positive class.
        control: The control label for the negative class.
            When there are only two classes in `in.metafile` at column `envs.clscol`,
            either `case` or `control` can be specified and the other will be
            automatically set to the other class.
        gmtfile: The pathways in GMT format, with the gene names/ids in the same format as the seurat object.
            One could also use a URL to a GMT file. For example, from <https://download.baderlab.org/EM_Genesets/current_release/Human/symbol/Pathways/>.
        method (choice): The method to do the preranking.
            - signal_to_noise: Signal to noise.
                The larger the differences of the means (scaled by the standard deviations);
                that is, the more distinct the gene expression is in each phenotype and the more the gene
                acts as a "class marker".
            - s2n: Alias of signal_to_noise.
            - abs_signal_to_noise: The absolute value of signal_to_noise.
            - abs_s2n: Alias of abs_signal_to_noise.
            - t_test: T test.
                Uses the difference of means scaled by the standard deviation and number of samples.
            - ratio_of_classes: Also referred to as fold change.
                Uses the ratio of class means to calculate fold change for natural scale data.
            - diff_of_classes: Difference of class means.
                Uses the difference of class means to calculate fold change for nature scale data
            - log2_ratio_of_classes: Log2 ratio of class means.
                Uses the log2 ratio of class means to calculate fold change for natural scale data.
                This is the recommended statistic for calculating fold change for log scale data.
        clscol: The column of metafile specifying the classes of the samples
            When `in.metafile` is not specified, it can also be specified as a list of
            classes, in the same order as the samples in `in.infile`.
        top (type=auto): Do gsea table and enrich plot for top N pathways.
            If it is < 1, will apply it to `padj`, selecting pathways with `padj` < `top`.
        eps (type=float): This parameter sets the boundary for calculating the p value.
            See <https://rdrr.io/bioc/fgsea/man/fgseaMultilevel.html>
        minsize (type=int): Minimal size of a gene set to test. All pathways below the threshold are excluded.
        maxsize (type=int): Maximal size of a gene set to test. All pathways above the threshold are excluded.
        rest (type=json;order=98): Rest arguments for [`fgsea()`](https://rdrr.io/bioc/fgsea/man/fgsea.html)
            See also <https://rdrr.io/bioc/fgsea/man/fgseaMultilevel.html>
        cases (type=json;order=99): If you have multiple cases, you can specify them here.
            The keys are the names of the cases and the values are the above options except `mutaters`.
            If some options are not specified, the default values specified above will be used.
            If no cases are specified, the default case will be added with the name `GSEA`.

    Requires:
        bioconductor-fgsea:
            - check: {{proc.lang}} -e "library(fgsea)"
    """  # noqa: E501
    input = "infile:file, metafile:file"
    output = "outdir:dir:{{in.infile | stem}}.fgsea"
    lang = config.lang.rscript
    envs = {
        "ncores": config.misc.ncores,
        "case": None,
        "control": None,
        "gmtfile": None,
        "method": "signal_to_noise",
        "clscol": None,
        "top": 10,
        "eps": 0,
        "minsize": 10,
        "maxsize": 100,
        "rest": {},
        "cases": {},
    }
    script = "file://../scripts/gsea/FGSEA.R"
    plugin_opts = {"report": "file://../reports/gsea/FGSEA.svelte"}


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
    plugin_opts = {"report": "file://../reports/gsea/Enrichr.svelte"}
