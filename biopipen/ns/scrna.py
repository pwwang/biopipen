"""Tools to analyze single-cell RNA"""

from ..core.proc import Proc
from ..core.config import config


class SeuratLoading(Proc):
    """Seurat - Loading data

    Deprecated, should be superseded by SeuratPreparing

    Input:
        metafile: The metadata of the samples
            A tab-delimited file
            Two columns are required:
            - `Sample` to specify the sample names.
            - `RNADir` to assign the path of the data to the samples
              The path will be read by `Read10X()` from `Seurat`

    Output:
        rdsfile: The RDS file with a list of Seurat object

    Envs:
        qc: The QC filter for each sample.
            This will be passed to `subset(obj, subset=<qc>)`.
            For example
            `nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5`
    """

    input = "metafile:file"
    output = "rdsfile:file:{{in.metafile | stem}}.seurat.RDS"
    envs = {"qc": ""}
    lang = config.lang.rscript
    script = "file://../scripts/scrna/SeuratLoading.R"


class SeuratPreparing(Proc):
    """Seurat - Loading, preparing and applying QC to data

    This process will -
    - Prepare the seurat object
    - Apply QC to the data

    See also -
    https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#standard-pre-processing-workflow-1)
    and
    https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Create_one_merged_object

    Input:
        metafile: The metadata of the samples
            A tab-delimited file
            Two columns are required:
            `Sample` to specify the sample names.
            `RNADir` to assign the path of the data to the samples
            The path will be read by `Read10X()` from `Seurat`

    Output:
        rdsfile: The RDS file with the Seurat object
            Note that the cell ids are preficed with sample names
            QC plots will be saved in `<job.outdir>/before-qc` and
            `<job.outdir>/after-qc`

    Envs:
        ncores: Number of cores to use
        cell_qc: filter expression to filter cells, using
            `tidyrseurat::filter()`.
            Available QC keys include `nFeature_RNA`, `nCount_RNA`,
            `percent.mt`, `percent.ribo`, `percent.hb`, and `percent.plat`
            For example: `nFeature_RNA > 200 & percent.mt < 5` will
            keep cells with more than 200 genes and less than 5%% mitochondrial
            genes.
        gene_qc (ns): filter genes. Currently only `min_cells` is supported.
            `gene_qc` is applied after `cell_qc`.
            - min_cells: the minimum number of cells that a gene must be
                expressed in to be kept.

    Requires:
        r-seurat:
            - check: {{proc.lang}} <(echo "library(Seurat)")
        r-future:
            - check: {{proc.lang}} <(echo "library(future)")
        r-bracer:
            - check: {{proc.lang}} <(echo "library(bracer)")
    """  # noqa: E501
    input = "metafile:file"
    output = "rdsfile:file:{{in.metafile | stem}}.seurat.RDS"
    lang = config.lang.rscript
    envs = {
        "ncores": config.misc.ncores,
        "cell_qc": None,  # "nFeature_RNA > 200 & percent.mt < 5",
        "gene_qc": {"min_cells": 3},
    }
    script = "file://../scripts/scrna/SeuratPreparing.R"
    plugin_opts = {
        "report": "file://../reports/scrna/SeuratPreparing.svelte",
    }


class SeuratClustering(Proc):
    """Seurat - Determine the clusters

    Input:
        srtobj: The seurat object loaded by SeuratPreparing

    Output:
        rdsfile: The seurat object with cluster information

    Envs:
        ncores: Number of cores to use
        use_sct (action=store_true): Whether use SCTransform routine or not
            See https://satijalab.org/seurat/articles/integration_rpca.html#performing-integration-on-datasets-normalized-with-sctransform-1
            and https://satijalab.org/seurat/articles/integration_rpca.html
        SCTransform (ns): Arguments for `SCTransform()`
            - method: The method to use for normalization.
            - <more>: See https://satijalab.org/seurat/reference/sctransform
        SelectIntegrationFeatures (ns): Arguments for `SelectIntegrationFeatures()`
            - nfeatures: The number of features to select
            - <more>: See https://satijalab.org/seurat/reference/selectintegrationfeatures
        PrepSCTIntegration (ns): Arguments for `PrepSCTIntegration()`
            - <more>: See https://satijalab.org/seurat/reference/prepsctintegration
        NormalizeData (ns): Arguments for `NormalizeData()`
            - <more>: See https://satijalab.org/seurat/reference/normalizedata
        FindVariableFeatures (ns): Arguments for `FindVariableFeatures()`
            - <more>: See https://satijalab.org/seurat/reference/findvariablefeatures
        FindIntegrationAnchors (ns): Arguments for `FindIntegrationAnchors()`
            - <more>: See https://satijalab.org/seurat/reference/findintegrationanchors
        IntegrateData (ns): Arguments for `IntegrateData()`
            - <more>: See https://satijalab.org/seurat/reference/integratedata
        ScaleData (ns): Arguments for `ScaleData()`
            - verbose (action=store_true): Whether to print the progress
            - <more>: See https://satijalab.org/seurat/reference/scaledata
        RunPCA (ns): Arguments for `RunPCA()`
            - npcs (type=int): The number of PCs to compute
            - verbose (action=store_true): Whether to print the progress
            - <more>: See https://satijalab.org/seurat/reference/runpca
        RunUMAP (ns): Arguments for `RunUMAP()`
            - dims (type=int): The number of PCs to use
            - reduction: The reduction to use for UMAP
            - <more>: See https://satijalab.org/seurat/reference/runumap
        FindNeighbors (ns): Arguments for `FindNeighbors()`
            - <more>: See https://satijalab.org/seurat/reference/findneighbors
        FindClusters (ns): Arguments for `FindClusters()`
            - resolution: The resolution of the clustering
            - <more>: See https://satijalab.org/seurat/reference/findclusters

    Requires:
        r-seurat:
            - check: {{proc.lang}} <(echo "library(Seurat)")
        r-tidyr:
            - check: {{proc.lang}} <(echo "library(tidyr)")
        r-dplyr:
            - check: {{proc.lang}} <(echo "library(dplyr)")
    """  # noqa: E501
    input = "srtobj:file"
    output = "rdsfile:file:{{in.srtobj | stem}}.RDS"
    lang = config.lang.rscript
    envs = {
        "ncores": config.misc.ncores,
        "use_sct": False,
        "SCTransform": {"method": "glmGamPoi"},
        "SelectIntegrationFeatures": {"nfeatures": 3000},
        "PrepSCTIntegration": {},
        "NormalizeData": {},
        "FindVariableFeatures": {},
        "FindIntegrationAnchors": {},
        "IntegrateData": {},
        "ScaleData": {"verbose": False},
        "RunPCA": {"npcs": 30, "verbose": False},
        "RunUMAP": {"reduction": "pca", "dims": 30},
        "FindNeighbors": {},
        "FindClusters": {"resolution": 0.8},
    }
    script = "file://../scripts/scrna/SeuratClustering.R"


class SeuratClusterStats(Proc):
    """Seurat - Cluster statistics

    Input:
        srtobj: The seurat object loaded by SeuratClustering

    Output:
        outdir: The output directory

    Envs:
        stats (ctype=json): The statistics to plot
            `nCells_*` - Number of cells for each cluster.
            You can specify `by` to group the cells by a metadata column,
            and `devpars` to specify the device parameters for the plot.
            You can also specify `filter` to filter the cells under certain
            conditions using metadata columns.
            `fracCells_*` - Fraction of cells for each cluster.
            Similar to `nCells_*`, but the fraction is calculated
            instead of the number.
        exprs (ctype=json): The expression values to plot.
            `genes` - The set of genes for the plots, unless `features` for
            those plots is specified. Could also specify a file with genes
            (one per line).
            `ridgeplots` - The ridge plots for the gene expressions.
            See `?Seurat::RidgePlot`.
            `vlnplots` - Violin plots for the gene expressions.
            See `?Seurat::VlnPlot`. You can have `boxplot` key to add
            `geom_boxplot()` to the violin plots.
            `featureplots` - The feature plots for the gene expressions.
            See `?Seurat::FeaturePlot`.
            `dotplot` - Dot plots for the gene expressions.
            See `?Seurat::DotPlot`.
            `heatmap` - Heatmap for the gene expressions.
            See `?Seurat::DoHeatmap`. You can specify `average=True` to plot on
            the average of the expressions.
            `table` - The table for the gene expressions
            (supported keys: title, log2, subset and features).
            All the above with `devpars` to define the output figures
            and `plus` to add elements to the `ggplot` object.
            You can have `subset` to subset the data. Multiple cases can be
            distinguished by `ridgeplots` and `ridgeplots_1`
            If no `features` specified, will use `genes`. If you want to use
            the default gene list `VariantFeatures(srtobj)[1:20]`, specify
            `features = "default"`. Or you can also specify the genes directly
            to `features`.
        dimplots (ctype=json): The dimensional reduction plots.
            `<case>` - The case to plot. Keys are the arguments for
            `Seurat::Dimplot()`, add `devpars`.

    Requires:
        r-seurat:
            - check: {{proc.lang}} -e "library(Seurat)"
    """

    input = "srtobj:file"
    output = "outdir:dir:{{in.srtobj | stem}}.cluster_stats"
    lang = config.lang.rscript
    envs = {
        "stats": {
            "ncells_All": {},
            "ncells_Sample": {"by": "Sample"},
            "fracCells_Sample": {"by": "Sample"},
        },
        "exprs": {},
        "dimplots": {
            "Ident": {
                "group-by": "ident",
                "devpars": {"res": 100, "height": 1000, "width": 1000},
            }
        },
    }
    script = "file://../scripts/scrna/SeuratClusterStats.R"
    plugin_opts = {
        "report": "file://../reports/scrna/SeuratClusterStats.svelte"
    }


class CellsDistribution(Proc):
    """Distribution of cells (i.e. in a TCR clone) from different groups
    for each cluster

    This generates a set of pie charts with proportion of cells in each cluster
    Rows are the cells identities (i.e. TCR clones or TCR clusters), columns
    are groups (i.e. clinic groups).

    Input:
        srtobj: The seurat object generated by SeuratClustering
        casefile: The file with the cases
            See `envs.cases`

    Output:
        outdir: The output directory

    Envs:
        name: The name of the job.
        cases (ctype=json): The cases to use.
            If `in.casefile` is not provided, `envs.name` and `envs.cases`
            will be used.
            >>> # The name of the job, used in report, optional
            >>> # If not given, will use `{{in.srtobj | stem}}`, ...
            >>> name = ""
            >>> [cases.case1]
            >>> filter = 'Responder != "Other"'
            >>> # Create some helper columns
            >>> [cases.case1.mutaters]
            >>> Responder = '''
            >>>     case_when(
            >>>         PFS == "CONTROL" ~ "CONTROL",
            >>>         PFS >= 60 ~ "Responder",
            >>>         PFS < 60 ~ "Non-Responder",
            >>>         TRUE ~ "Other"
            >>>     )
            >>> '''
            >>> CloneSize = '''
            >>>     case_when(
            >>>         Clones > 50 ~ "Clone50p",
            >>>         Clones > 30 & Clones <= 50 ~ "Clone30p",
            >>>         Clones > 10 & Clones <= 30 ~ "Clone10p",
            >>>         TRUE ~ "CloneSmall"
            >>>     )
            >>> '''
            >>> [cases.case1.group]
            >>> # Use values of a column from metadata for the groups
            >>> by = "Timepoint"
            >>> # Specify the order of the groups to show on the plot
            >>> # from left to right (columns)
            >>> # Optional
            >>> order = ["CONTROL", "Responder", "Non-Responder"]
            >>> #
            >>> [cases.case1.cells]
            >>> by = "CDR3.aa" # "CloneSize"
            >>> # Specify the order of the cell identities to show on the plot
            >>> # The order of the cells (rows)
            >>> orderby = "desc(CloneSize)"
            >>> n = 10


    Requires:
        r-seurat:
            - check: {{proc.lang}} -e "library(Seurat)"
        r-dplyr:
            - check: {{proc.lang}} -e "library(dplyr)"
        r-tidyr:
            - check: {{proc.lang}} -e "library(tidyr)"
    """
    input = "srtobj:file, casefile:file"
    output = "outdir:dir:{{in.srtobj | stem}}.cells_distribution"
    lang = config.lang.rscript
    envs = {"name": None, "cases": {}}
    script = "file://../scripts/scrna/CellsDistribution.R"
    plugin_opts = {"report": "file://../reports/scrna/CellsDistribution.svelte"}


class SeuratMetadataMutater(Proc):
    """Mutate the metadata of the seurat object

    Input:
        srtobj: The seurat object loaded by SeuratPreparing
        metafile: Additional metadata
            A tab-delimited file with columns as meta columns and rows as
            cells.
        mutaters: A TOML string, file or a python dict to create new columns
            in metadata. They key-value will be evaluated
            with `srtobj@meta.data |> mutate(<key> = <value>)`
            The keys could be from `in.metafile`

    Output:
        rdsfile: The seurat object with the additional metadata

    Requires:
        r-seurat:
            - check: {{proc.lang}} <(echo "library(Seurat)")
        r-tibble:
            - check: {{proc.lang}} <(echo "library(tibble)")
        r-dplyr:
            - check: {{proc.lang}} <(echo "library(dplyr)")
    """
    input = "srtobj:file, metafile:file, mutaters:var"
    output = "rdsfile:file:{{in.srtobj | stem}}.RDS"
    lang = config.lang.rscript
    script = "file://../scripts/scrna/SeuratMetadataMutater.R"


class GeneExpressionInvestigation(Proc):
    """Investigation of expressions of genes of interest

    Input:
        srtobj: The seurat object loaded by `SeuratPreparing`
        genefile: The genes to show their expressions in the plots
            Either one column or two columns.
            If one column, the column name will be used as both the gene names
            to match the expressions and the names to show in the plots
            If two columns, the first column will be used as the gene names
            to match the expressions and the second column will be used to
            show in the plots.
        configfile: The configuration file (toml). See `envs.config`
            If not provided, use `envs.config`

    Output:
        outdir: The output directory with the plots

    Envs:
        gopts: Options for `read.table()` to read `in.genefile`
        config: The configurations to do the plots
            name: The name of the job, mostly used in report
            mutaters: The mutater to mutate the metadata
            groupby: Which meta columns to group the data
            subset: Select a subset of cells, will be passed to
                `subset(obj, subset=<subset>)`
            plots: Plots to generate
                Currently supported
                `boxplot`:
                - `ncol`: Split the plot to how many columns?
                - `res`, `height` and `width` the parameters for `png()`
                `heatmap`:
                - `res`, `height` and `width` the parameters for `png()`
                - other arguments for `ComplexHeatmap::Heatmap()`
    """
    input = "srtobj:file, genefile:file, configfile:file"
    output = "outdir:dir:{{in.configfile | stem0}}.gei"
    lang = config.lang.rscript
    order = 4
    envs = {
        "config": {},
        "gopts": {
            "header": False,
            "row.names": None,
            "sep": "\t",
            "check.names": False,
        },
    }
    script = "file://../scripts/scrna/GeneExpressionInvistigation.R"
    plugin_opts = {
        "report": "file://../reports/scrna/GeneExpressionInvistigation.svelte"
    }


class DimPlots(Proc):
    """Seurat - Dimensional reduction plots

    Input:
        srtobj: The seruat object in RDS format
        configfile: A toml configuration file with "cases"
            If this is given, `envs.cases` will be overriden
        name: The name of the job, used in report

    Output:
        outdir: The output directory

    Envs:
        cases: The cases for the dim plots
            Keys are the names and values are the arguments to
            `Seurat::Dimplots`
    """
    input = "srtobj:file, configfile:file, name:var"
    output = "outdir:dir:{{in.srtobj | stem}}.dimplots"
    lang = config.lang.rscript
    script = "file://../scripts/scrna/DimPlots.R"
    envs = {"cases": {"Ident": {"group.by": "ident"}}}
    plugin_opts = {
        "report": "file://../reports/scrna/DimPlots.svelte",
        "report_toc": False,
    }


class MarkersFinder(Proc):
    """Find markers between different groups of cells

    Input:
        srtobj: The seurat object loaded by `SeuratPreparing`
        casefile: The config file in TOML that looks like. See `envs.cases`

    Output:
        outdir: The output directory for the markers

    Envs:
        ncores: Number of cores to use to parallelize the groups
        cases (ctype=json): The cases to find markers for.
            >>> # The name of the job, used in report
            >>> name = ""
            >>> [cases.case1]
            >>> "ident-1" = "Tumor"
            >>> "ident-2" = "Normal"
            >>> "group-by" = "Source"
            >>> # focus on a subset of cells
            >>> filter = "SampleType != 'Control'"
            >>> # other arguments for Seruat::FindMarkers()
            We can also use a new group.by:
            >>> [cases.case2]
            >>> "ident-1" = "Case"
            >>> "ident-2" = "Control"
            >>> "group-by" = "Group"
            >>> # Do the comparison in each cluster
            >>> each = "seurat_clusters"
            >>> # other arguments for Seruat::FindMarkers()
            >>> # Filter after mutaters
            >>> filter2 = "SampleType != 'Control'"
            >>> [cases.case2.mutaters]
            >>> Group = '''
            >>>   if_else(Source %%in%% c("Tumor", "Normal"), "Case", "Control")
            >>> '''
            If "ident-2" is not provided, it will use the rest of the cells
            as "ident-2".
            If only "group.by" is given, will call `FindAllMarkers()`
        dbs (ctype=list): The dbs to do enrichment analysis for significant
            markers See below for all librarys
            https://maayanlab.cloud/Enrichr/#libraries
        sigmarkers: An expression passed to `dplyr::filter()` to filter the
            significant markers for enrichment analysis
            Available variables are `p_val`, `avg_log2FC`, `pct.1`, `pct.2` and
            `p_val_adj`
            Example - `"p_val_adj < 0.05 & abs(avg_log2FC) > 1"`
    """
    input = "srtobj:file, casefile:file"
    output = "outdir:dir:{{(in.casefile or in.srtobj) | stem0}}.markers"
    lang = config.lang.rscript
    envs = {
        "ncores": config.misc.ncores,
        "cases": {
            "Cluster": {"group-by": "seurat_clusters"},
        },
        "dbs": [
            "GO_Biological_Process_2021",
            "GO_Cellular_Component_2021",
            "GO_Molecular_Function_2021",
            "KEGG_2021_Human",
        ],
        "sigmarkers": "p_val_adj < 0.05",
    }
    order = 5
    script = "file://../scripts/scrna/MarkersFinder.R"
    plugin_opts = {"report": "file://../reports/scrna/MarkersFinder.svelte"}


class ExprImpute(Proc):
    """Impute the dropout values in scRNA-seq data.

    Input:
        infile: The input file in RDS format of Seurat object

    Output:
        outfile: The output file in RDS format of Seurat object
            Note that with rmagic and alra, the original RNA assay will be
            renamed to `UNIMPUTED_RNA` and the imputed RNA assay will be
            renamed to `RNA`

    Envs:
        tool (choice): Either alra, scimpute or rmagic
            - alra: Use RunALRA() from Seurat
            - scimpute: Use scImpute() from scimpute
            - rmagic: Use magic() from Rmagic
        scimpute_args (ns): The arguments for scimpute
            - drop_thre (type=float): The dropout threshold
            - kcluster (type=int): Number of clusters to use
            - ncores (type=int): Number of cores to use
            - refgene: The reference gene file
        rmagic_args (ns): The arguments for rmagic
            - python: The python path where magic-impute is installed.
        alra_args (ctype=json): The arguments for `RunALRA()`

    Requires:
        r-scimpute:
            - if: {{proc.envs.tool == "scimpute"}}
            - check: {{proc.lang}} <(echo "library(scImpute)")
        r-rmagic:
            - if: {{proc.envs.tool == "rmagic"}}
            - check: |
                {{proc.lang}} <(\
                    echo "\
                        tryCatch(\
                            { setwd(dirname(Sys.getenv('CONDA_PREFIX'))) }, \
                            error = function(e) NULL \
                        ); \
                        library(Rmagic)\
                    "\
                )
        magic-impute:
            - if: {{proc.envs.tool == "rmagic"}}
            - check: {{proc.envs.rmagic_args.python}} -c "import magic")
        r-dplyr:
            - if: {{proc.envs.tool == "scimpute"}}
            - check: {{proc.lang}} <(echo "library(dplyr)")
        r-seurat:
            - check: {{proc.lang}} <(echo "library(Seurat)")
        r-seuratwrappers:
            - if: {{proc.envs.tool == "alra"}}
            - check: {{proc.lang}} <(echo "library(SeuratWrappers)")
    """
    input = "infile:file"
    output = "outfile:file:{{in.infile | stem}}.imputed.RDS"
    lang = config.lang.rscript
    envs = {
        "tool": "alra",
        "rmagic_args": {"python": config.exe.magic_python},
        "scimpute_args": {
            "drop_thre": 0.5,
            "kcluster": None,
            "ncores": config.misc.ncores,
            "refgene": config.ref.refgene,
        },
        "alra_args": {},
    }
    script = "file://../scripts/scrna/ExprImpute.R"


class SCImpute(Proc):
    """Impute the dropout values in scRNA-seq data.

    Deprecated. Use `ExprImpute` instead.

    Input:
        infile: The input file for imputation
            Either a SeuratObject or a matrix of count/TPM
        groupfile: The file to subset the matrix or label the cells
            Could be an output from ImmunarchFilter

    Output:
        outfile: The output matrix

    Envs:
        infmt: The input format.
            Either `seurat` or `matrix`
    """
    input = "infile:file, groupfile:file"
    output = [
        "outfile:file:{{in.infile | stem | replace: '.seurat', ''}}."
        "{{envs.outfmt}}"
    ]
    lang = config.lang.rscript
    envs = {
        "infmt": "seurat",  # or matrix
        "outfmt": "txt",  # or csv, rds
        "type": "count",  # or TPM
        "drop_thre": 0.5,
        "kcluster": None,  # or Kcluster
        "ncores": config.misc.ncores,
        "refgene": config.ref.refgene,
    }
    script = "file://../scripts/scrna/SCImpute.R"


class SeuratFilter(Proc):
    """Filtering cells from a seurat object

    Input:
        srtobj: The seurat object in RDS
        filters: The filters to apply. Could be a file or string in TOML, or
            a python dictionary, with following keys:
            - mutaters: Create new columns in the metadata
            - filter: A R expression that will pass to
              `subset(sobj, subset = ...)` to filter the cells

    Output:
        outfile: The filtered seurat object in RDS

    Envs:
        invert: Invert the selection?

    Requires:
        r-seurat:
            - check: {{proc.lang}} <(echo "library('Seurat')")
        r-dplyr:
            - check: {{proc.lang}} <(echo "library('dplyr')")
    """
    input = "srtobj:file, filters:var"
    output = "outfile:file:{{in.srtobj | stem}}.filtered.RDS"
    lang = config.lang.rscript
    envs = {"invert": False}
    script = "file://../scripts/scrna/SeuratFilter.R"


class SeuratSubset(Proc):
    """Subset a seurat object into multiple seruat objects

    Input:
        srtobj: The seurat object in RDS
        subsets: The subsettings to apply. Could be a file or string in TOML, or
            a python dictionary, with following keys:
            - <name>: Name of the case
                mutaters: Create new columns in the metadata
                subset: A R expression that will pass to
                    `subset(sobj, subset = ...)`
                groupby: The column to group by, each value will be a case
                    If groupby is given, subset will be ignored, each value
                    of the groupby column will be a case

    Output:
        outdir: The output directory with the subset seurat objects

    Envs:
        ignore_nas: Ignore NA values?

    Requires:
        r-seurat:
            - check: {{proc.lang}} <(echo "library('Seurat')")
        r-dplyr:
            - check: {{proc.lang}} <(echo "library('dplyr')")
    """
    input = "srtobj:file, subsets:var"
    output = "outdir:dir:{{in.srtobj | stem}}.subsets"
    envs = {"ignore_nas": True}
    lang = config.lang.rscript
    script = "file://../scripts/scrna/SeuratSubset.R"


class SeuratSplit(Proc):
    """Split a seurat object into multiple seruat objects

    Input:
        srtobj: The seurat object in RDS
        by: The metadata column to split by

    Output:
        outdir: The output directory with the subset seurat objects

    Envs:
        by: The metadata column to split by
            Ignored if `by` is given in the input
        recell: Rename the cell ids using the `by` column
            A string of R function taking the original cell ids and `by`
    """
    input = "srtobj:file, by:var"
    output = "outdir:dir:{{in.srtobj | stem}}.subsets"
    envs = {
        "by": None,
        "recell": None,  # "function(cellid, by) {}",
    }
    lang = config.lang.rscript
    script = "file://../scripts/scrna/SeuratSplit.R"


class Subset10X(Proc):
    """Subset 10X data, mostly used for testing

    Requires r-matrix to load matrix.mtx.gz

    Input:
        indir: The input directory

    Output:
        outdir: The output directory

    Envs:
        seed: The seed for random number generator
        nfeats: The number of features to keep.
            If <=1 then it will be the percentage of features to keep
        ncells: The number of cells to keep.
            If <=1 then it will be the percentage of cells to keep
        feats_to_keep: The features/genes to keep.
            The final features list will be `feats_to_keep` + `nfeats`
    """
    input = "indir:dir"
    output = "outdir:dir:{{in.indir | stem}}"
    envs = {
        "seed": 8525,
        "nfeats": 0.1,
        "ncells": 0.1,
        "feats_to_keep": [],
    }
    lang = config.lang.rscript
    script = "file://../scripts/scrna/Subset10X.R"


class Write10X(Proc):
    """Write a Seurat object to 10X format

    using `write10xCounts` from `DropletUtils`

    Input:
        srtobj: The seurat object in RDS

    Output:
        outdir: The output directory

    Envs:
        version: The version of 10X format
    """
    input = "srtobj:file"
    output = "outdir:dir:{{in.srtobj | stem}}"
    envs = {"version": "3"}
    lang = config.lang.rscript
    script = "file://../scripts/scrna/Write10X.R"


class ScFGSEA(Proc):
    """Fast gene set enrichment analysis (fgsea) for cells in different groups

    Input:
        srtobj: The seurat object loaded by `SeuratPreparing`
        casefile: The config file in TOML
            See `MarkersFinder.envs.cases`

    Output:
        outdir: The output directory for the results

    Envs:
        name: The name of the job, used in report
        ncores: Number of cores to use to parallelize the groups
        cases (ctype=json): The cases to find markers for.
            `ident-2` is required for each case.
            One could also use placeholders for the cases.
            To enable this, you need `percluster = True` in the config.
            Currently only cluster is supported. One could use `{cluster}` or
            `{ident}` to denote the clusters.
        gmtfile: The pathways in GMT format
        method (choice): The method to do the preranking.
            Signal to noise: the larger the differences of the means
            (scaled by the standard deviations); that is, the more distinct
            the gene expression is in each phenotype and the more the gene
            acts as a “class marker.”.
            Absolute signal to noise: the absolute value of the signal to noise.
            T test: Uses the difference of means scaled by the standard
            deviation and number of samples.
            Ratio of classes: Uses the ratio of class means to calculate
            fold change for natural scale data.
            Diff of classes: Uses the difference of class means to calculate
            fold change for nature scale data
            Log2 ratio of classes: Uses the log2 ratio of class means to
            calculate fold change for natural scale data. This is the
            recommended statistic for calculating fold change for log scale
            data.
            - signal_to_noise: Signal to noise
            - s2n: Alias of signal_to_noise
            - abs_signal_to_noise: absolute signal to noise
            - abs_s2n: Alias of abs_signal_to_noise
            - t_test: T test
            - ratio_of_classes: Also referred to as fold change
            - diff_of_classes: Difference of class means
            - log2_ratio_of_classes: Log2 ratio of class means
        top (ctype=auto): Do gsea table and enrich plot for top N pathways.
            If it is < 1, will apply it to `padj`
        <rest>: Rest arguments for `fgsea()`

    Requires:
        bioconductor-fgsea:
            - check: {{proc.lang}} -e "library(fgsea)"
        r-seurat:
            - check: {{proc.lang}} -e "library(seurat)"
    """
    input = "srtobj:file, casefile:file"
    output = "outdir:dir:{{(in.casefile or in.srtobj) | stem0}}.fgsea"
    lang = config.lang.rscript
    envs = {
        "name": None,
        "ncores": config.misc.ncores,
        "cases": {},
        "gmtfile": "",
        "method": "s2n",
        "top": 20,
        "minSize": 10,
        "maxSize": 100,
        "eps": 0,
    }
    script = "file://../scripts/scrna/ScFGSEA.R"
    plugin_opts = {"report": "file://../reports/scrna/ScFGSEA.svelte"}


class CellTypeAnnotate(Proc):
    """Annotate cell types

    Input:
        sobjfile: The seurat object

    Output:
        outfile: The rds file of seurat object with cell type annotated

    Envs:
        tool (choice): The tool to use for cell type annotation.
            sctype: https://github.com/IanevskiAleksandr/sc-type .
            scCATCH: https://github.com/ZJUFanLab/scCATCH
            - sctype: Use `scType` to annotate cell types
            - sccatch: Use `scCATCH` to annotate cell types
            - direct: Directly assign cell types
        sctype_tissue: The tissue to use for sctype.
            E.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,
            Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus
        sctype_db: The database to use for sctype
        cell_types (list): The cell types to use for direct annotation
            Each a list of cell type names, or a dict with keys as the old
            identity and values as the new cell type.

    Requires:
        r-HGNChelper:
            - if: {{proc.envs.tool == 'sctype'}}
            - check: {{proc.lang}} -e "library(HGNChelper)"
        r-seurat:
            - if: {{proc.envs.tool == 'sctype'}}
            - check: {{proc.lang}} -e "library(Seurat)"
        r-dplyr:
            - if: {{proc.envs.tool == 'sctype'}}
            - check: {{proc.lang}} -e "library(dplyr)"
        r-openxlsx:
            - if: {{proc.envs.tool == 'sctype'}}
            - check: {{proc.lang}} -e "library(openxlsx)"
    """
    input = "sobjfile:file"
    output = "outfile:file:{{in.sobjfile | stem}}.annotated.RDS"
    lang = config.lang.rscript
    envs = {
        "tool": "sctype",
        "sctype_tissue": None,
        "sctype_db": config.ref.sctype_db,
        "cell_types": [],
    }
    script = "file://../scripts/scrna/CellTypeAnnotate.R"


class SeuratMap2Ref(Proc):
    """Map the seurat object to reference

    See: https://satijalab.org/seurat/articles/integration_mapping.html
    and https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

    Input:
        sobjfile: The seurat object

    Output:
        outfile: The rds file of seurat object with cell type annotated

    Envs:
        use (choice): Which level of cell type to use for further analysis and
            being aliased to `alias`
            - predicted.celltype.l1: The first level of predicted cell type
            - predicted.celltype.l2: The second level of predicted cell type
        alias: The name of an aliasied column to `use`
            This is helpful for the downstream analysis where the column name
            is used as the cluster.
        ref: The reference seurat object file.
            Either an RDS file or a h5seurat file that can be loaded by
            `Seurat::LoadH5Seurat()`.
            The file type is determined by the extension. `.rds` or `.RDS` for
            RDS file, `.h5seurat` or `.h5` for h5seurat file.
        SCTransform (ns): Arguments for `SCTransform()`
            - <more>: See https://satijalab.org/seurat/reference/sctransform
        FindTransferAnchors (ns): Arguments for `FindTransferAnchors()`
            - <more>: See https://satijalab.org/seurat/reference/findtransferanchors
        MapQuery (ns): Arguments for `MapQuery()`
            - refdata (ctype=json): Data to transfer
            - <more>: See https://satijalab.org/seurat/reference/mapquery

    Requires:
        r-seurat:
            - check: {{proc.lang}} -e "library(Seurat)"
    """  # noqa: E501
    input = "sobjfile:file"
    output = "outfile:file:{{in.sobjfile | stem}}.RDS"
    lang = config.lang.rscript
    envs = {
        "use": "predicted.celltype.l2",
        "alias": "seurat_clusters",
        "ref": None,
        "SCTransform": {
            "do-correct-umi": False,
            "do-scale": False,
            "do-center": True,
        },
        "FindTransferAnchors": {
            "normalization-method": "SCT",
            "reference-reduction": "spca",
        },
        "MapQuery": {
            "reference-reduction": "spca",
            "reduction-model": "wnn.umap",
            "refdata": {
                "celltype-l1": "celltype.l1",
                "celltype-l2": "celltype.l2",
                "predicted_ADT": "ADT",
            }
        },
    }
    script = "file://../scripts/scrna/SeuratMap2Ref.R"
