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
    """Seurat - Loading and preparing data

    What will be done ?
    https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#standard-pre-processing-workflow-1)
    1. All samples with be integrated as a single seurat object
    2. QC
    3. Normalization
    4. Feature selection
    5. Scaling
    6. Linear dimensional reduction

    Input:
        metafile: The metadata of the samples
            A tab-delimited file
            Two columns are required:
            - `Sample` to specify the sample names.
            - `RNADir` to assign the path of the data to the samples
              The path will be read by `Read10X()` from `Seurat`

    Output:
        rdsfile: The RDS file with the Seurat object
            Note that the cell ids are preficed with sample names

    Envs:
        ncores: Number of cores to use

    """

    input = "metafile:file"
    output = "rdsfile:file:{{in.metafile | stem}}.seurat.RDS"
    lang = config.lang.rscript
    envs = {
        "ncores": config.misc.ncores,
    }
    script = "file://../scripts/scrna/SeuratPreparing.R"
    # plugin_opts = {
    #     "report": "file://../reports/scrna/SeuratPreparing.svelte"
    # }


class SeuratClustering(Proc):
    """Seurat - Determine the clusters

    Input:
        srtobj: The seurat object loaded by SeuratPreparing

    Output:
        rdsfile: The seurat object with cluster information
        groupfile: A groupfile with cells for future analysis

    Envs:
        FindClusters: Arguments to `FindClusters()`
    """

    input = "srtobj:file"
    output = [
        "rdsfile:file:{{in.srtobj | stem}}.RDS",
        "groupfile:file:{{in.srtobj | stem | replace: '.seurat', ''}}"
        ".clusters.txt",
    ]
    lang = config.lang.rscript
    envs = {"FindClusters": {"resolution": 0.8}}
    script = "file://../scripts/scrna/SeuratClustering.R"


class GeneExpressionInvestigation(Proc):
    """Investigation of expressions of genes of interest

    Input:
        srtobj: The seurat object loaded by SeuratPreparing
        groupfile: The group of cells with the first column the groups and
            rest the cells in each sample.
            Or the subset conditions using metadata of `srtobj`
            See `envs.group_subset`
        genefiles: The genes to show their expressions in the plots
        configfile: The configuration file (toml). See `envs`
            If not provided, use `envs`
    Output:
        outdir: The output directory with the plots

    Envs:
        group_subset: Is the `in.groupfile` subset conditions using metadata
            Or the groupfile as described.
        name: The name to name the job. Otherwise the stem of groupfile
            will be used
        target` Which sample to pull expression from could be multiple
        gopts: Options for `read.table()` to read the genefiles
        plots: Plots to generate for this case
            `boxplot`:
            - `use`: Which gene file to use (1-based)
            - `ncol`: Split the plot to how many columns?
            - `res`, `height` and `width` the parameters for `png()`
            `heatmap`:
            - `use`: Which gene file to use (1-based)
            - `res`, `height` and `width` the parameters for `png()`
            - other arguments for `ComplexHeatmap::Heatmap()`
    """

    input = "srtobj:file, groupfile:file, genefiles:files, configfile:file"
    output = "outdir:dir:{{in.configfile | stem0}}.exprs"
    lang = config.lang.rscript
    order = 4
    envs = {
        "group_subset": False,
        "name": None,
        "target": None,
        "gopts": {},
        "plots": {},
    }
    script = "file://../scripts/scrna/GeneExpressionInvistigation.R"
    plugin_opts = {
        "report": "file://../reports/scrna/GeneExpressionInvistigation.svelte"
    }


class DimPlots(Proc):
    """Seurat - Dimensional reduction plots"""

    input = "srtobj:file, groupfile:file, configfile:file"
    output = "outfile:file:{{in.groupfile | stem}}.dimplot.png"
    lang = config.lang.rscript
    script = "file://../scripts/scrna/DimPlots.R"
    plugin_opts = {
        "report": "file://../reports/scrna/DimPlots.svelte",
        "report_toc": False,
    }


class MarkersFinder(Proc):
    """Find markers between different groups of cells

    Input:
        srtobj: The seurat object loaded by SeuratLoading
        groupfile: The group of cells with first column the groups and
            rest the cells in each sample.
        casefile: The config file in TOML that looks like
            >>> [case1]
            >>> "ident.1" = "ident.1"
            >>> "ident.2" = "ident.2"
            >>> # other arguments for Seruat::FindMarkers()
        name: The name of the jobs, mosted used in report

    Output:
        outdir: The output directory for the markers

    Envs:
        ncores: Number of cores to use to parallelize the groups
        cases: The cases to find markers for.
            Values would be the arguments for `FindMarkers()`
            If "ALL" or "ALL" in the keys, the process will run for all groups
            in the groupfile. The other keys will be arguments to `FindMarkers`
            When `ident.2` is not given and there is only one group or more
            than two groups in groupfile, the rest cells in the object will
            be used as the control
            if it is `ident`, will igore the groupfile and find markers for all
            idents.
        dbs: The dbs to do enrichment analysis for significant markers
            See below for all librarys
            https://maayanlab.cloud/Enrichr/#libraries
    """

    input = "srtobj:file, groupfile:file, name:var, casefile:file"
    output = "outdir:dir:{{in.groupfile | stem0}}.markers"
    lang = config.lang.rscript
    envs = {
        "ncores": config.misc.ncores,
        "cases": {"ALL": True},
        "dbs": [
            "GO_Biological_Process_2021",
            "GO_Cellular_Component_2021",
            "GO_Molecular_Function_2021",
            "KEGG_2021_Human",
        ],
    }
    order = 5
    script = "file://../scripts/scrna/MarkersFinder.R"
    plugin_opts = {"report": "file://../reports/scrna/MarkersFinder.svelte"}


class SCImpute(Proc):
    """Impute the dropout values in scRNA-seq data.

    Input:
        infile: The input file for imputation
            Either a SeuratObject or a matrix of count/TPM
        groupfile: The file to subset the matrix or label the cells
            Could be an output from ImmunarchFilter

    Output:
        outfile: The output matrix

    Envs:
        infmt: The input format.
            Either `seurat` or `matrix
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
        filterfile: The file with the filtering information
            Either a group file (rows cases for filtering, columns are samples
            or ALL for all cells with prefices), or config under
            `subsetting` section in TOML with keys
            and value that will be passed to `subset(...,subset = ...)`

    Output:
        out: The filtered seurat object in RDS if `envs.multicase` is False,
            otherwise the directory with the filtered seurat objects

    Envs:
        filterfmt: `auto`, `subset` or `grouping`.
            If `subset` then `in.filterfile` will be config in TOML, otherwise
            if `grouping`, it is a groupfile. See `in.filterfile`.
            If `auto`, test if there is `=` in the file. If so, it's `subset`
            otherwise `grouping`
        invert: Invert the selection?
        multicase: If True, multiple seurat objects will be generated.
            For `envs.filterfmt == "subset"`, each key-value pair will be a
            case, otherwise, each row of `in.filterfile` will be a case.

    """
    input = "srtobj:file, filterfile:file"
    output = """
        {%- if envs.multicase -%}
            out:dir:{{in.filterfile | stem}}.seuratfiltered
        {%- elif envs.filterfmt == "subset" -%}
            out:file:{{in.filterfile | read | toml_loads | list | first}}.RDS
        {%- else -%}
            out:file:{{in.filterfile | readlines | last | split | first}}.RDS
        {%- endif -%}
    """
    envs = {
        "filterfmt": "auto",  # subset or grouping
        "invert": False,
        "multicase": True,
    }
    lang = config.lang.rscript
    script = "file://../scripts/scrna/SeuratFilter.R"
