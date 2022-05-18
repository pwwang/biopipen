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

    Envs:
        FindClusters: Arguments to `FindClusters()`
    """

    input = "srtobj:file"
    output = "rdsfile:file:{{in.srtobj | stem}}.RDS"
    lang = config.lang.rscript
    envs = {"FindClusters": {"resolution": 0.8}}
    script = "file://../scripts/scrna/SeuratClustering.R"


class SeuratMetadataMutater(Proc):
    """Mutate the metadata of the seurat object

    Input:
        srtobj: The seurat object loaded by SeuratPreparing
        metafile: Additional metadata
            A tab-delimited file with columns as meta columns and rows as
            cells.

    Output:
        rdsfile: The seurat object with the additional metadata

    Envs:
        mutaters: A dictionary of mutates with keys the new meta columns
            and values the new meta values. They key-value will be evaluated
            with `srtobj@meta.data |> mutate(<key> = <value>)`
            The keys could be from `in.metafile`
    """
    input = "srtobj:file, metafile:file"
    output = "rdsfile:file:{{in.srtobj | stem}}.RDS"
    lang = config.lang.rscript
    envs = {"mutaters": {}}
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
    envs = {
        "cases": {
            "Ident": {"group.by": "ident"}
        }
    }
    plugin_opts = {
        "report": "file://../reports/scrna/DimPlots.svelte",
        "report_toc": False,
    }


class MarkersFinder(Proc):
    """Find markers between different groups of cells

    Input:
        srtobj: The seurat object loaded by `SeuratPreparing`
        casefile: The config file in TOML that looks like

            >>> # The name of the job, used in report
            >>> name = ""
            >>> [cases.case1]
            >>> "ident.1" = "Tumor"
            >>> "ident.2" = "Normal"
            >>> "group.by" = "Source"
            >>> # other arguments for Seruat::FindMarkers()

            We can also use a new group.by:
            >>> [cases.case2]
            >>> "ident.1" = "Case"
            >>> "ident.2" = "Control"
            >>> "group.by" = "Group"
            >>> # other arguments for Seruat::FindMarkers()
            >>> [cases.case2.mutaters]
            >>> Group = '''
            >>>   if_else(Source %in% c("Tumor", "Normal"), "Case", "Control")
            >>> '''

            If "ident.2" is not provided, it will use the rest of the cells
            as "ident.2".

            If only "group.by" is given, will call `FindAllMarkers()`

    Output:
        outdir: The output directory for the markers

    Envs:
        ncores: Number of cores to use to parallelize the groups
        cases: The cases to find markers for.
            See `in.casefile`.
        dbs: The dbs to do enrichment analysis for significant markers
            See below for all librarys
            https://maayanlab.cloud/Enrichr/#libraries
    """

    input = "srtobj:file, casefile:file"
    output = "outdir:dir:{{(in.casefile or in.srtobj) | stem0}}.markers"
    lang = config.lang.rscript
    envs = {
        "ncores": config.misc.ncores,
        "cases": {
            "name": "Markers for all clusters",
            "cases": {
                "Cluster": {"group.by": "seurat_clusters"},
            }
        },
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
