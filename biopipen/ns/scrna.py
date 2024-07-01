"""Tools to analyze single-cell RNA"""

from ..core.proc import Proc
from ..core.config import config
from ..utils.common_docstrs import (
    indent_docstr,
    format_placeholder,
    MUTATE_HELPERS_CLONESIZE,
    ENVS_SECTION_EACH,
)

MUTATE_HELPERS_CLONESIZE_INDENTED = indent_docstr(MUTATE_HELPERS_CLONESIZE, "    " * 3)
ENVS_SECTION_EACH_INDENTED = indent_docstr(ENVS_SECTION_EACH, "    " * 3)


class SeuratLoading(Proc):
    """Seurat - Loading data

    Deprecated, should be superseded by SeuratPreparing

    Input:
        metafile: The metadata of the samples
            A tab-delimited file
            Two columns are required:
            - `Sample` to specify the sample names.
            - `RNAData` to assign the path of the data to the samples
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
    """Load, prepare and apply QC to data, using `Seurat`

    This process will -
    - Prepare the seurat object
    - Apply QC to the data
    - Integrate the data from different samples

    See also
    - <https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#standard-pre-processing-workflow-1)>
    - <https://satijalab.org/seurat/articles/integration_introduction>

    This process will read the scRNA-seq data, based on the information provided by
    `SampleInfo`, specifically, the paths specified by the `RNAData` column.
    Those paths should be either paths to directoies containing `matrix.mtx`,
    `barcodes.tsv` and `features.tsv` files that can be loaded by
    [`Seurat::Read10X()`](https://satijalab.org/seurat/reference/read10x),
    or paths to `h5` files that can be loaded by
    [`Seurat::Read10X_h5()`](https://satijalab.org/seurat/reference/read10x_h5).

    Each sample will be loaded individually and then merged into one `Seurat` object, and then perform QC.

    In order to perform QC, some additional columns are added to the meta data of the `Seurat` object. They are:

    - `precent.mt`: The percentage of mitochondrial genes.
    - `percent.ribo`: The percentage of ribosomal genes.
    - `precent.hb`: The percentage of hemoglobin genes.
    - `percent.plat`: The percentage of platelet genes.

    For integration, two routes are available:

    - [Performing integration on datasets normalized with `SCTransform`](https://satijalab.org/seurat/articles/seurat5_integration#perform-streamlined-one-line-integrative-analysis)
    - [Using `NormalizeData` and `FindIntegrationAnchors`](https://satijalab.org/seurat/articles/seurat5_integration#layers-in-the-seurat-v5-object)

    /// Note
    When using `SCTransform`, the default Assay will be set to `SCT` in output, rather than `RNA`.
    If you are using `cca` or `rpca` interation, the default assay will be `integrated`.
    ///

    /// Note
    From `biopipen` v0.23.0, this requires `Seurat` v5.0.0 or higher.
    ///

    Input:
        metafile: The metadata of the samples
            A tab-delimited file
            Two columns are required:
            `Sample` to specify the sample names.
            `RNAData` to assign the path of the data to the samples
            The path will be read by `Read10X()` from `Seurat`, or the path
            to the h5 file that can be read by `Read10X_h5()` from `Seurat`.

    Output:
        rdsfile: The RDS file with the Seurat object with all samples integrated.
            Note that the cell ids are preficed with sample names QC plots will be
            saved in `<job.outdir>/before-qc` and `<job.outdir>/after-qc`.

    Envs:
        ncores (type=int): Number of cores to use.
            Used in `future::plan(strategy = "multicore", workers = <ncores>)`
            to parallelize some Seurat procedures.
        cell_qc: Filter expression to filter cells, using
            `tidyrseurat::filter()`.
            Available QC keys include `nFeature_RNA`, `nCount_RNA`,
            `percent.mt`, `percent.ribo`, `percent.hb`, and `percent.plat`.

            /// Tip | Example
            Including the columns added above, all available QC keys include
            `nFeature_RNA`, `nCount_RNA`, `percent.mt`, `percent.ribo`, `percent.hb`,
            and `percent.plat`. For example:

            ```toml
            [SeuratPreparing.envs]
            cell_qc = "nFeature_RNA > 200 & percent.mt < 5"
            ```
            will keep cells with more than 200 genes and less than 5%% mitochondrial
            genes.
            ///

        cell_qc_per_sample (flag): Whether to perform cell QC per sample or not.
            If `True`, the cell QC will be performed per sample, and the QC will be
            applied to each sample before merging.
        gene_qc (ns): Filter genes.
            `gene_qc` is applied after `cell_qc`.
            - min_cells: The minimum number of cells that a gene must be
                expressed in to be kept.
            - excludes: The genes to exclude. Multiple genes can be specified by
                comma separated values, or as a list.

            /// Tip | Example
            ```toml
            [SeuratPreparing.envs]
            gene_qc = { min_cells = 3 }
            ```
            will keep genes that are expressed in at least 3 cells.
            ///

        use_sct (flag): Whether use SCTransform routine to integrate samples or not.
            Before the following procedures, the `RNA` layer will be split by samples.

            If `False`, following procedures will be performed in the order:
            * [`NormalizeData`](https://satijalab.org/seurat/reference/normalizedata).
            * [`FindVariableFeatures`](https://satijalab.org/seurat/reference/findvariablefeatures).
            * [`ScaleData`](https://satijalab.org/seurat/reference/scaledata).
            See <https://satijalab.org/seurat/articles/seurat5_integration#layers-in-the-seurat-v5-object>
            and <https://satijalab.org/seurat/articles/pbmc3k_tutorial.html>

            If `True`, following procedures will be performed in the order:
            * [`SCTransform`](https://satijalab.org/seurat/reference/sctransform).
            See <https://satijalab.org/seurat/articles/seurat5_integration#perform-streamlined-one-line-integrative-analysis>

        no_integration (flag): Whether to skip integration or not.
        NormalizeData (ns): Arguments for [`NormalizeData()`](https://satijalab.org/seurat/reference/normalizedata).
            `object` is specified internally, and `-` in the key will be replaced with `.`.
            - <more>: See <https://satijalab.org/seurat/reference/normalizedata>

        FindVariableFeatures (ns): Arguments for [`FindVariableFeatures()`](https://satijalab.org/seurat/reference/findvariablefeatures).
            `object` is specified internally, and `-` in the key will be replaced with `.`.
            - <more>: See <https://satijalab.org/seurat/reference/findvariablefeatures>

        ScaleData (ns): Arguments for [`ScaleData()`](https://satijalab.org/seurat/reference/scaledata).
            `object` and `features` is specified internally, and `-` in the key will be replaced with `.`.
            - <more>: See <https://satijalab.org/seurat/reference/scaledata>

        RunPCA (ns): Arguments for [`RunPCA()`](https://satijalab.org/seurat/reference/runpca).
            `object` and `features` is specified internally, and `-` in the key will be replaced with `.`.
            - npcs (type=int): The number of PCs to compute.
                For each sample, `npcs` will be no larger than the number of columns - 1.
            - <more>: See <https://satijalab.org/seurat/reference/runpca>

        SCTransform (ns): Arguments for [`SCTransform()`](https://satijalab.org/seurat/reference/sctransform).
            `object` is specified internally, and `-` in the key will be replaced with `.`.
            - `return-only-var-genes`: Whether to return only variable genes.
            - `min_cells`: The minimum number of cells that a gene must be expressed in to be kept.
                A hidden argument of `SCTransform` to filter genes.
                If you try to keep all genes in the `RNA` assay, you can set `min_cells` to `0` and
                `return-only-var-genes` to `False`.
                See <https://github.com/satijalab/seurat/issues/3598#issuecomment-715505537>
            - <more>: See <https://satijalab.org/seurat/reference/sctransform>

        IntegrateLayers (ns): Arguments for [`IntegrateLayers()`](https://satijalab.org/seurat/reference/integratelayers).
            `object` is specified internally, and `-` in the key will be replaced with `.`.
            When `use_sct` is `True`, `normalization-method` defaults to `SCT`.
            - method (choice): The method to use for integration.
                - CCAIntegration: Use `Seurat::CCAIntegration`.
                - CCA: Same as `CCAIntegration`.
                - cca: Same as `CCAIntegration`.
                - RPCAIntegration: Use `Seurat::RPCAIntegration`.
                - RPCA: Same as `RPCAIntegration`.
                - rpca: Same as `RPCAIntegration`.
                - HarmonyIntegration: Use `Seurat::HarmonyIntegration`.
                - Harmony: Same as `HarmonyIntegration`.
                - harmony: Same as `HarmonyIntegration`.
                - FastMNNIntegration: Use `Seurat::FastMNNIntegration`.
                - FastMNN: Same as `FastMNNIntegration`.
                - fastmnn: Same as `FastMNNIntegration`.
                - scVIIntegration: Use `Seurat::scVIIntegration`.
                - scVI: Same as `scVIIntegration`.
                - scvi: Same as `scVIIntegration`.
            - <more>: See <https://satijalab.org/seurat/reference/integratelayers>

        DoubletFinder (ns): Arguments to run [`DoubletFinder`](https://github.com/chris-mcginnis-ucsf/DoubletFinder).
            See also <https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/DoubletFinder.html>.
            To disable `DoubletFinder`, set `envs.DoubletFinder` to `None` or `False`; or set `pcs` to `0`.
            - PCs (type=int): Number of PCs to use for 'doubletFinder' function.
            - doublets (type=float): Number of expected doublets as a proportion of the pool size.
            - pN (type=float): Number of doublets to simulate as a proportion of the pool size.
            - ncores (type=int): Number of cores to use for `DoubletFinder::paramSweep`.
                Set to `None` to use `envs.ncores`.
                Since parallelization of the function usually exhausts memory, if big `envs.ncores` does not work
                for `DoubletFinder`, set this to a smaller number.

        cache (type=auto): Whether to cache the information at different steps.
            If `True`, the seurat object will be cached in the job output directory, which will be not cleaned up when job is rerunning.
            The cached seurat object will be saved as `<signature>.<kind>.RDS` file, where `<signature>` is the signature determined by
            the input and envs of the process.
            See <https://github.com/satijalab/seurat/issues/7849>, <https://github.com/satijalab/seurat/issues/5358> and
            <https://github.com/satijalab/seurat/issues/6748> for more details also about reproducibility issues.
            To not use the cached seurat object, you can either set `cache` to `False` or delete the cached file at
            `<signature>.RDS` in the cache directory.

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
        "cell_qc_per_sample": False,
        "gene_qc": {"min_cells": 0, "excludes": []},
        "use_sct": False,
        "no_integration": False,
        "NormalizeData": {},
        "FindVariableFeatures": {},
        "ScaleData": {},
        "RunPCA": {},
        "SCTransform": {
            "return-only-var-genes": True,
            "min_cells": 5,
        },
        "IntegrateLayers": {"method": "harmony"},
        "DoubletFinder": {"PCs": 0, "pN": 0.25, "doublets": 0.075, "ncores": 1},
        "cache": config.path.tmpdir,
    }
    script = "file://../scripts/scrna/SeuratPreparing.R"
    plugin_opts = {
        "report": "file://../reports/scrna/SeuratPreparing.svelte",
    }


class SeuratClustering(Proc):
    """Determine the clusters of cells without reference using Seurat FindClusters
    procedure.

    Input:
        srtobj: The seurat object loaded by SeuratPreparing

    Output:
        rdsfile: The seurat object with cluster information at `seurat_clusters`
            If `SCTransform` was used, the default Assay will be reset to `RNA`.

    Envs:
        ncores (type=int;order=-100): Number of cores to use.
            Used in `future::plan(strategy = "multicore", workers = <ncores>)`
            to parallelize some Seurat procedures.
            See also: <https://satijalab.org/seurat/articles/future_vignette.html>
        ScaleData (ns): Arguments for [`ScaleData()`](https://satijalab.org/seurat/reference/scaledata).
            If you want to re-scale the data by regressing to some variables, `Seurat::ScaleData`
            will be called. If nothing is specified, `Seurat::ScaleData` will not be called.
            - vars-to-regress: The variables to regress on.
            - <more>: See <https://satijalab.org/seurat/reference/scaledata>
        SCTransform (ns): Arguments for [`SCTransform()`](https://satijalab.org/seurat/reference/sctransform).
            If you want to re-scale the data by regressing to some variables, `Seurat::SCTransform`
            will be called. If nothing is specified, `Seurat::SCTransform` will not be called.
            - vars-to-regress: The variables to regress on.
            - <more>: See <https://satijalab.org/seurat/reference/sctransform>
        RunUMAP (ns): Arguments for [`RunUMAP()`](https://satijalab.org/seurat/reference/runumap).
            `object` is specified internally, and `-` in the key will be replaced with `.`.
            `dims=N` will be expanded to `dims=1:N`; The maximal value of `N` will be the minimum of `N` and the number of columns - 1 for each sample.
            - dims (type=int): The number of PCs to use
            - reduction: The reduction to use for UMAP.
                If not provided, `sobj@misc$integrated_new_reduction` will be used.
            - <more>: See <https://satijalab.org/seurat/reference/runumap>
        FindNeighbors (ns): Arguments for [`FindNeighbors()`](https://satijalab.org/seurat/reference/findneighbors).
            `object` is specified internally, and `-` in the key will be replaced with `.`.
            - reduction: The reduction to use.
                If not provided, `sobj@misc$integrated_new_reduction` will be used.
            - <more>: See <https://satijalab.org/seurat/reference/findneighbors>
        FindClusters (ns): Arguments for [`FindClusters()`](https://satijalab.org/seurat/reference/findclusters).
            `object` is specified internally, and `-` in the key will be replaced with `.`.
            The cluster labels will be saved in `seurat_clusters` and prefixed with "c".
            The first cluster will be "c1", instead of "c0".
            - resolution (type=auto): The resolution of the clustering. You can have multiple resolutions as a list or as a string separated by comma.
                Ranges are also supported, for example: `0.1:0.5:0.1` will generate `0.1, 0.2, 0.3, 0.4, 0.5`. The step can be omitted, defaulting to 0.1.
                The results will be saved in `seurat_clusters_<resolution>`.
                The final resolution will be used to define the clusters at `seurat_clusters`.
            - <more>: See <https://satijalab.org/seurat/reference/findclusters>
        clustree_devpars (ns): The device parameters for the clustree plots.
            - res (type=int): The resolution of the plots.
            - height (type=int): The height of the plots.
            - width (type=int): The width of the plots.
        cache (type=auto): Whether to cache the information at different steps.
            If `True`, the seurat object will be cached in the job output directory, which will be not cleaned up when job is rerunning.
            The cached seurat object will be saved as `<signature>.<kind>.RDS` file, where `<signature>` is the signature determined by
            the input and envs of the process.
            See <https://github.com/satijalab/seurat/issues/7849>, <https://github.com/satijalab/seurat/issues/5358> and
            <https://github.com/satijalab/seurat/issues/6748> for more details also about reproducibility issues.
            To not use the cached seurat object, you can either set `cache` to `False` or delete the cached file at
            `<signature>.RDS` in the cache directory.

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
        "ScaleData": {},
        "SCTransform": {},
        "RunUMAP": {"dims": 30},
        "FindNeighbors": {},
        "FindClusters": {"resolution": 0.8},
        "clustree_devpars": {"res": 100, "height": 1000, "width": 800},
        "cache": config.path.tmpdir,
    }
    script = "file://../scripts/scrna/SeuratClustering.R"


class SeuratSubClustering(Proc):
    """Find clusters of a subset of cells.

    It's unlike [`Seurat::FindSubCluster`], which only finds subclusters of a single
    cluster. Instead, it will perform the whole clustering procedure on the subset of
    cells. One can use metadata to specify the subset of cells to perform clustering on.

    For the subset of cells, the reductions will be re-performed on the subset of cells,
    and then the clustering will be performed on the subset of cells. The reduction
    will be saved in `sobj@reduction$sub_umap_<casename>` of the original object and the
    clustering will be saved in the metadata of the original object using the casename \
    as the column name.

    Input:
        srtobj: The seurat object

    Output:
        rdsfile: The seurat object with the subclustering information.

    Envs:
        ncores (type=int;order=-100): Number of cores to use.
            Used in `future::plan(strategy = "multicore", workers = <ncores>)`
            to parallelize some Seurat procedures.
        mutaters (type=json): The mutaters to mutate the metadata to subset the cells.
            The mutaters will be applied in the order specified.
        subset: An expression to subset the cells, will be passed to
            [`tidyseurat::filter()`](https://stemangiola.github.io/tidyseurat/reference/filter.html).

        RunUMAP (ns): Arguments for [`RunUMAP()`](https://satijalab.org/seurat/reference/runumap).
            `object` is specified internally as the subset object, and `-` in the key will be replaced with `.`.
            `dims=N` will be expanded to `dims=1:N`; The maximal value of `N` will be the minimum of `N` and the number of columns - 1 for each sample.
            - dims (type=int): The number of PCs to use
            - reduction: The reduction to use for UMAP.
                If not provided, `sobj@misc$integrated_new_reduction` will be used.
            - <more>: See <https://satijalab.org/seurat/reference/runumap>
        FindNeighbors (ns): Arguments for [`FindNeighbors()`](https://satijalab.org/seurat/reference/findneighbors).
            `object` is specified internally, and `-` in the key will be replaced with `.`.
            - reduction: The reduction to use.
                If not provided, `sobj@misc$integrated_new_reduction` will be used.
            - <more>: See <https://satijalab.org/seurat/reference/findneighbors>
        FindClusters (ns): Arguments for [`FindClusters()`](https://satijalab.org/seurat/reference/findclusters).
            `object` is specified internally, and `-` in the key will be replaced with `.`.
            The cluster labels will be prefixed with "s". The first cluster will be "s1", instead of "s0".
            - resolution (type=auto): The resolution of the clustering. You can have multiple resolutions as a list or as a string separated by comma.
                Ranges are also supported, for example: `0.1:0.5:0.1` will generate `0.1, 0.2, 0.3, 0.4, 0.5`. The step can be omitted, defaulting to 0.1.
                The results will be saved in `<casename>_<resolution>`.
                The final resolution will be used to define the clusters at `<casename>`.
            - <more>: See <https://satijalab.org/seurat/reference/findclusters>
        clustree_devpars (ns): The device parameters for the clustree plots.
            - res (type=int): The resolution of the plots.
            - height (type=int): The height of the plots.
            - width (type=int): The width of the plots.
        cache (type=auto): Whether to cache the information at different steps.
            If `True`, the seurat object will be cached in the job output directory, which will be not cleaned up when job is rerunning.
            The cached seurat object will be saved as `<signature>.<kind>.RDS` file, where `<signature>` is the signature determined by
            the input and envs of the process.
            See <https://github.com/satijalab/seurat/issues/7849>, <https://github.com/satijalab/seurat/issues/5358> and
            <https://github.com/satijalab/seurat/issues/6748> for more details also about reproducibility issues.
            To not use the cached seurat object, you can either set `cache` to `False` or delete the cached file at
            `<signature>.RDS` in the cache directory.
        cases (type=json): The cases to perform subclustering.
            Keys are the names of the cases and values are the dicts inherited from `envs` except `mutaters` and `cache`.
            If empty, a case with name `subcluster` will be created with default parameters.
    """  # noqa: E501
    input = "srtobj:file"
    output = "rdsfile:file:{{in.srtobj | stem}}.RDS"
    lang = config.lang.rscript
    envs_depth = 1
    envs = {
        "ncores": config.misc.ncores,
        "mutaters": {},
        "subset": None,
        "RunUMAP": {"dims": 30},
        "FindNeighbors": {},
        "FindClusters": {"resolution": 0.8},
        "clustree_devpars": {"res": 100, "height": 1000, "width": 800},
        "cache": config.path.tmpdir,
        "cases": {"subcluster": {}},
    }
    script = "file://../scripts/scrna/SeuratSubClustering.R"


class SeuratClusterStats(Proc):
    """Statistics of the clustering.

    Including the number/fraction of cells in each cluster, the gene expression values
    and dimension reduction plots. It's also possible to perform stats on
    TCR clones/clusters or other metadata for each T-cell cluster.

    Examples:
        ### Number of cells in each cluster

        ```toml
        [SeuratClusterStats.envs.stats]
        # suppose you have nothing set in `envs.stats_defaults`
        # otherwise, the settings will be inherited here
        nCells_All = { }
        ```

        ![nCells_All](https://pwwang.github.io/immunopipe/latest/processes/images/SeuratClusterStats_nCells_All.png){: width="80%" }

        ### Number of cells in each cluster by groups

        ```toml
        [SeuratClusterStats.envs.stats]
        nCells_Sample = { group-by = "Sample" }
        ```

        ![nCells_Sample](https://pwwang.github.io/immunopipe/latest/processes/images/SeuratClusterStats_nCells_Sample.png){: width="80%" }

        ### Violin plots for the gene expressions

        ```toml
        [SeuratClusterStats.envs.features]
        features = "CD4,CD8A"
        # Remove the dots in the violin plots
        vlnplots = { pt-size = 0, kind = "vln" }
        # Don't use the default genes
        vlnplots_1 = { features = ["FOXP3", "IL2RA"], pt-size = 0, kind = "vln" }
        ```

        ![vlnplots](https://pwwang.github.io/immunopipe/latest/processes/images/SeuratClusterStats_vlnplots.png){: width="80%" }
        ![vlnplots_1](https://pwwang.github.io/immunopipe/latest/processes/images/SeuratClusterStats_vlnplots_1.png){: width="80%" }

        ### Dimension reduction plot with labels

        ```toml
        [SeuratClusterStats.envs.dimplots.Idents]
        label = true
        label-box = true
        repel = true
        ```

        ![dimplots](https://pwwang.github.io/immunopipe/latest/processes/images/SeuratClusterStats_dimplots.png){: width="80%" }

    Input:
        srtobj: The seurat object loaded by `SeuratClustering`

    Output:
        outdir: The output directory

    Envs:
        mutaters (type=json): The mutaters to mutate the metadata to subset the cells.
            The mutaters will be applied in the order specified.
        hists_defaults (ns): The default parameters for histograms.
            This will plot histograms for the number of cells along `x`.
            For example, you can plot the number of cells along cell activity score.
            - x: The column name in metadata to plot as the x-axis.
                The NA values will be removed.
                It could be either numeric or factor/character.
            - x_order (list): The order of the x-axis, only works for factor/character `x`.
                You can also use it to subset `x` (showing only a subset values of `x`).
            - cells_by: A column name in metadata to group the cells.
                The NA values will be removed. It should be a factor/character.
                if not specified, all cells will be used.
            - cells_order (list): The order of the cell groups for the plots.
                It should be a list of strings. You can also use `cells_orderby` and `cells_n`
                to determine the order.
            - cells_orderby: An expression passed to `dplyr::arrange()` to order the cell groups.
            - cells_n: The number of cell groups to show.
                Ignored if `cells_order` is specified.
            - ncol (type=int): The number of columns for the plots, split by `cells_by`.
            - subset: An expression to subset the cells, will be passed to `dplyr::filter()`.
            - each: Whether to plot each group separately.
            - bins: The number of bins to use, only works for numeric `x`.
            - plus (list): The extra elements to add to the `ggplot` object.
            - devpars (ns): The device parameters for the plots.
                - res (type=int): The resolution of the plots.
                - height (type=int): The height of the plots.
                - width (type=int): The width of the plots.
        hists (type=json): The cases for histograms.
            Keys are the names of the plots and values are the dicts inherited from `env.hists_defaults`.
            There is no default case.
        stats_defaults (ns): The default parameters for `stats`.
            This is to do some basic statistics on the clusters. For more comprehensive analysis,
            see `RadarPlots` and `CellsDistribution`.
            The parameters from the cases can overwrite the default parameters.
            - frac (flag): Whether to output the fraction of cells instead of number.
            - pie (flag): Also output a pie chart?
            - circos (flag): Also output a circos plot?
            - table (flag): Whether to output a table (in tab-delimited format) and in the report.
            - frac_ofall (flag): Whether to output the fraction against all cells,
                instead of the fraction in each group.
                Does not work for circos plot.
                Only works when `frac` is `True` and `group-by` is specified.
            - transpose (flag): Whether to transpose the cluster and group, that is,
                using group as the x-axis and cluster to fill the plot.
                For circos plot, when transposed, the arrows will be drawn from the idents (by `ident`) to the
                the groups (by `group-by`).
                Only works when `group-by` is specified.
            - position (choice): The position of the bars. Does not work for pie and circos plots.
                - stack: Use `position_stack()`.
                - fill: Use `position_fill()`.
                - dodge: Use `position_dodge()`.
                - auto: Use `stack` when there are more than 5 groups, otherwise use `dodge`.
            - ident: The column name in metadata to use as the identity.
            - group-by: The column name in metadata to group the cells.
                Does NOT support for pie charts.
            - split-by: The column name in metadata to split the cells into different plots.
                Does NOT support for circos plots.
            - subset: An expression to subset the cells, will be passed to
                `dplyr::filter()` on metadata.
            - circos_labels_rot (flag): Whether to rotate the labels in the circos plot.
                In case the labels are too long.
            - circos_devpars (ns): The device parameters for the circos plots.
                - res (type=int): The resolution of the plots.
                - height (type=int): The height of the plots.
                - width (type=int): The width of the plots.
            - pie_devpars (ns): The device parameters for the pie charts.
                - res (type=int): The resolution of the plots.
                - height (type=int): The height of the plots.
                - width (type=int): The width of the plots.
            - devpars (ns): The device parameters for the plots.
                - res (type=int): The resolution of the plots.
                - height (type=int): The height of the plots.
                - width (type=int): The width of the plots.
        stats (type=json): The number/fraction of cells to plot.
            Keys are the names of the plots and values are the dicts inherited from `env.stats_defaults`.
            Here are some examples -
            >>> {
            >>>     "nCells_All": {},
            >>>     "nCells_Sample": {"group-by": "Sample"},
            >>>     "fracCells_Sample": {"frac": True, "group-by": "Sample"},
            >>> }
        ngenes_defaults (ns): The default parameters for `ngenes`.
            The default parameters to plot the number of genes expressed in each cell.
            - ident: The column name in metadata to use as the identity.
            - group-by: The column name in metadata to group the cells.
                Dodge position will be used to separate the groups.
            - split-by: The column name in metadata to split the cells into different plots.
            - subset: An expression to subset the cells, will be passed to `tidyrseurat::filter()`.
            - devpars (ns): The device parameters for the plots.
                - res (type=int): The resolution of the plots.
                - height (type=int): The height of the plots.
                - width (type=int): The width of the plots.
        ngenes (type=json): The number of genes expressed in each cell.
            Keys are the names of the plots and values are the dicts inherited from `env.ngenes_defaults`.
        features_defaults (ns): The default parameters for `features`.
            - features: The features to plot.
                It can be either a string with comma separated features, a list of features, a file path with `file://` prefix with features
                (one per line), or an integer to use the top N features from `VariantFeatures(srtobj)`.
            - ident: The column name in metadata to use as the identity.
                If it is from subclustering (reduction `sub_umap_<ident>` exists), the reduction will be used.
            - cluster_orderby (type=auto): The order of the clusters to show on the plot.
                An expression passed to `dplyr::summarise()` on the grouped data frame (by `seurat_clusters`).
                The summary stat will be passed to `dplyr::arrange()` to order the clusters. It's applied on the whole meta.data before grouping and subsetting.
                For example, you can order the clusters by the activation score of
                the cluster: `desc(mean(ActivationScore, na.rm = TRUE))`, suppose you have a column
                `ActivationScore` in the metadata.
                You may also specify the literal order of the clusters by a list of strings.
            - subset: An expression to subset the cells, will be passed to `tidyrseurat::filter()`.
            - devpars (ns): The device parameters for the plots. Does not work for `table`.
                - res (type=int): The resolution of the plots.
                - height (type=int): The height of the plots.
                - width (type=int): The width of the plots.
            - plus: The extra elements to add to the `ggplot` object. Does not work for `table`.
            - group-by: Group cells in different ways (for example, orig.ident). Works for `ridge`, `vln`, and `dot`.
                It also works for `feature` as `shape.by` being passed to [`Seurat::FeaturePlot`](https://satijalab.org/seurat/reference/featureplot).
            - split-by: The column name in metadata to split the cells into different plots.
                It works for `vln`, `feature`, and `dot`.
            - assay: The assay to use.
            - layer: The layer to use.
            - reduction: The reduction to use. Only works for `feature`.
            - section: The section to put the plot in the report.
                If not specified, the case title will be used.
            - ncol (type=int): The number of columns for the plots.
            - kind (choice): The kind of the plot or table.
                - ridge: Use `Seurat::RidgePlot`.
                - ridgeplot: Same as `ridge`.
                - vln: Use `Seurat::VlnPlot`.
                - vlnplot: Same as `vln`.
                - violin: Same as `vln`.
                - violinplot: Same as `vln`.
                - feature: Use `Seurat::FeaturePlot`.
                - featureplot: Same as `feature`.
                - dot: Use `Seurat::DotPlot`.
                - dotplot: Same as `dot`.
                - bar: Bar plot on an aggregated feature.
                    The features must be a single feature, which will be either an  existing feature or an expression
                    passed to `dplyr::summarise()` (grouped by `ident`) on the existing features to create a new feature.
                - barplot: Same as `bar`.
                - heatmap: Use `Seurat::DoHeatmap`.
                - avgheatmap: Plot the average expression of the features in each cluster as a heatmap.
                - table: The table for the features, only gene expressions are supported.
                    (supported keys: ident, subset, and features).
        features (type=json): The plots for features, include gene expressions, and columns from metadata.
            Keys are the titles of the cases and values are the dicts inherited from `env.features_defaults`. It can also have other parameters from
            each Seurat function used by `kind`. Note that for argument name with `.`, you should use `-` instead.
        dimplots_defaults (ns): The default parameters for `dimplots`.
            - ident: The identity to use.
                If it is from subclustering (reduction `sub_umap_<ident>` exists), this reduction will be used if `reduction`
                is set to `dim` or `auto`.
            - group-by: Same as `ident` if not specified, to define how the points are colored.
            - na_group: The group name for NA values, use `None` to ignore NA values.
            - split-by: The column name in metadata to split the cells into different plots.
            - shape-by: The column name in metadata to use as the shape.
            - subset: An expression to subset the cells, will be passed to `tidyrseurat::filter()`.
            - devpars (ns): The device parameters for the plots.
                - res (type=int): The resolution of the plots.
                - height (type=int): The height of the plots.
                - width (type=int): The width of the plots.
            - reduction (choice): Which dimensionality reduction to use.
                - dim: Use `Seurat::DimPlot`.
                    First searches for `umap`, then `tsne`, then `pca`.
                    If `ident` is from subclustering, `sub_umap_<ident>` will be used.
                - auto: Same as `dim`
                - umap: Use `Seurat::UMAPPlot`.
                - tsne: Use `Seurat::TSNEPlot`.
                - pca: Use `Seurat::PCAPlot`.
            - <more>: See <https://satijalab.org/seurat/reference/dimplot>
        dimplots (type=json): The dimensional reduction plots.
            Keys are the titles of the plots and values are the dicts inherited from `env.dimplots_defaults`. It can also have other parameters from
            [`Seurat::DimPlot`](https://satijalab.org/seurat/reference/dimplot).

    Requires:
        r-seurat:
            - check: {{proc.lang}} -e "library(Seurat)"
    """  # noqa: E501
    input = "srtobj:file"
    output = "outdir:dir:{{in.srtobj | stem}}.cluster_stats"
    lang = config.lang.rscript
    envs = {
        "mutaters": {},
        "hists_defaults": {
            "x": None,
            "x_order": [],
            "cells_by": None,
            "cells_order": [],
            "cells_orderby": None,
            "cells_n": 10,
            "subset": None,
            "ncol": 2,
            "each": None,
            "bins": 30,
            "plus": [],
            "devpars": {"res": 100, "height": None, "width": None},
        },
        "hists": {},
        "stats_defaults": {
            "frac": False,
            "pie": False,
            "circos": False,
            "table": False,
            "position": "auto",
            "frac_ofall": False,
            "transpose": False,
            "ident": "seurat_clusters",
            "group-by": None,
            "split-by": None,
            "subset": None,
            "circos_labels_rot": False,
            "devpars": {"res": 100, "height": 600, "width": 800},
            "pie_devpars": {"res": 100, "height": 600, "width": 800},
            "circos_devpars": {"res": 100, "height": 600, "width": 600},
        },
        "stats": {
            "Number of cells in each cluster": {
                "pie": True,
            },
            "Number of cells in each cluster by Sample": {
                "group-by": "Sample",
                "table": True,
                "frac": True,
            },
        },
        "ngenes_defaults": {
            "ident": "seurat_clusters",
            "group-by": None,
            "split-by": None,
            "subset": None,
            "devpars": {"res": 100, "height": 800, "width": 1000},
        },
        "ngenes": {
            "Number of genes expressed in each cluster": {},
        },
        "features_defaults": {
            "features": None,
            "ident": "seurat_clusters",
            "cluster_orderby": None,
            "subset": None,
            "devpars": {"res": 100},
            "plus": None,
            "group-by": None,
            "split-by": None,
            "assay": None,
            "section": None,
            "layer": None,
            "reduction": None,
            "kind": None,
            "ncol": 2,
        },
        "features": {},
        "dimplots_defaults": {
            "ident": "seurat_clusters",
            "group-by": None,
            "na_group": None,
            "split-by": None,
            "shape-by": None,
            "subset": None,
            "reduction": "dim",
            "devpars": {"res": 100, "height": 800, "width": 1000},
        },
        "dimplots": {
            "Dimensional reduction plot": {
                "label": True,
                "label-box": True,
                "repel": True,
            },
        },
    }
    script = "file://../scripts/scrna/SeuratClusterStats.R"
    plugin_opts = {
        "report": "file://../reports/scrna/SeuratClusterStats.svelte"
    }


class ModuleScoreCalculator(Proc):
    """Calculate the module scores for each cell

    The module scores are calculated by
    [`Seurat::AddModuleScore()`](https://satijalab.org/seurat/reference/addmodulescore)
    or [`Seurat::CellCycleScoring()`](https://satijalab.org/seurat/reference/cellcyclescoring)
    for cell cycle scores.

    The module scores are calculated as the average expression levels of each
    program on single cell level, subtracted by the aggregated expression of
    control feature sets. All analyzed features are binned based on averaged
    expression, and the control features are randomly selected from each bin.

    Input:
        srtobj: The seurat object loaded by `SeuratClustering`

    Output:
        rdsfile: The seurat object with module scores

    Envs:
        defaults (ns): The default parameters for `modules`.
            - features: The features to calculate the scores. Multiple features
                should be separated by comma.
                You can also specify `cc.genes` or `cc.genes.updated.2019` to
                use the cell cycle genes to calculate cell cycle scores.
                If so, three columns will be added to the metadata, including
                `S.Score`, `G2M.Score` and `Phase`.
                Only one type of cell cycle scores can be calculated at a time.
            - nbin (type=int): Number of bins of aggregate expression levels
                for all analyzed features.
            - ctrl (type=int): Number of control features selected from
                the same bin per analyzed feature.
            - k (flag): Use feature clusters returned from `DoKMeans`.
            - assay: The assay to use.
            - seed (type=int): Set a random seed.
            - search (flag): Search for symbol synonyms for features in
                features that don't match features in object?
            - keep (flag): Keep the scores for each feature?
                Only works for non-cell cycle scores.
            - agg (choice): The aggregation function to use.
                Only works for non-cell cycle scores.
                - mean: The mean of the expression levels
                - median: The median of the expression levels
                - sum: The sum of the expression levels
                - max: The max of the expression levels
                - min: The min of the expression levels
                - var: The variance of the expression levels
                - sd: The standard deviation of the expression levels
        modules (type=json): The modules to calculate the scores.
            Keys are the names of the expression programs and values are the
            dicts inherited from `env.defaults`.
            Here are some examples -
            >>> {
            >>>     "CellCycle": {"features": "cc.genes.updated.2019"},
            >>>     "Exhaustion": {"features": "HAVCR2,ENTPD1,LAYN,LAG3"},
            >>>     "Activation": {"features": "IFNG"},
            >>>     "Proliferation": {"features": "STMN1,TUBB"}
            >>> }

            For `CellCycle`, the columns `S.Score`, `G2M.Score` and `Phase` will
            be added to the metadata. `S.Score` and `G2M.Score` are the cell cycle
            scores for each cell, and `Phase` is the cell cycle phase for each cell.

            You can also add Diffusion Components (DC) to the modules
            >>> {"DC": {"features": 2, "kind": "diffmap"}}
            will perform diffusion map as a reduction and add the first 2
            components as `DC_1` and `DC_2` to the metadata. `diffmap` is a shortcut
            for `diffusion_map`. Other key-value pairs will pass to
            [`destiny::DiffusionMap()`](https://www.rdocumentation.org/packages/destiny/versions/2.0.4/topics/DiffusionMap%20class).
            You can later plot the diffusion map by using
            `reduction = "DC"` in `env.dimplots` in `SeuratClusterStats`.
            This requires [`SingleCellExperiment`](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
            and [`destiny`](https://bioconductor.org/packages/release/bioc/html/destiny.html) R packages.
    """  # noqa: E501
    input = "srtobj:file"
    output = "rdsfile:file:{{in.srtobj | stem}}.RDS"
    lang = config.lang.rscript
    envs = {
        "defaults": {
            "features": None,
            "nbin": 24,
            "ctrl": 100,
            "k": False,
            "assay": None,
            "seed": 8525,
            "search": False,
            "keep": False,
            "agg": "mean",
        },
        "modules": {
            # "CellCycle": {"features": "cc.genes.updated.2019"},
            # "Exhaustion": {"features": "HAVCR2,ENTPD1,LAYN,LAG3"},
            # "Activation": {"features": "IFNG"},
            # "Proliferation": {"features": "STMN1,TUBB"},
        },
    }
    script = "file://../scripts/scrna/ModuleScoreCalculator.R"


@format_placeholder(
    mutate_helpers_clonesize=MUTATE_HELPERS_CLONESIZE_INDENTED,
    envs_section_each=ENVS_SECTION_EACH_INDENTED,
)
class CellsDistribution(Proc):
    """Distribution of cells (i.e. in a TCR clone) from different groups
    for each cluster

    This generates a set of pie charts with proportion of cells in each cluster
    Rows are the cells identities (i.e. TCR clones or TCR clusters), columns
    are groups (i.e. clinic groups).

    Examples:
        ```toml
        [CellsDistribution.envs.mutaters]
        # Add Patient1_Tumor_Expanded column with CDR3.aa that
        # expands in Tumor of patient 1
        Patient1_Tumor_Expanded = '''
          expanded(., region, "Tumor", subset = patient == "Lung1", uniq = FALSE)
        '''

        [CellsDistribution.envs.cases.Patient1_Tumor_Expanded]
        cells_by = "Patient1_Tumor_Expanded"
        cells_orderby = "desc(CloneSize)"
        group_by = "region"
        group_order = [ "Tumor", "Normal" ]
        ```

        ![CellsDistribution_example](https://pwwang.github.io/immunopipe/latest/processes/images/CellsDistribution_example.png)

    Input:
        srtobj: The seurat object in RDS format

    Output:
        outdir: The output directory

    Envs:
        mutaters (type=json): The mutaters to mutate the metadata
            Keys are the names of the mutaters and values are the R expressions
            passed by `dplyr::mutate()` to mutate the metadata.
            %(mutate_helpers_clonesize)s

        cluster_orderby: The order of the clusters to show on the plot.
            An expression passed to `dplyr::summarise()` on the grouped data frame (by `seurat_clusters`).
            The summary stat will be passed to `dplyr::arrange()` to order the clusters. It's applied on the whole meta.data before grouping and subsetting.
            For example, you can order the clusters by the activation score of
            the cluster: `desc(mean(ActivationScore, na.rm = TRUE))`, suppose you have a column
            `ActivationScore` in the metadata.
        group_by: The column name in metadata to group the cells for the columns of the plot.
        group_order (list): The order of the groups (columns) to show on the plot
        cells_by: The column name in metadata to group the cells for the rows of the plot.
            If your cell groups have overlapping cells, you can also use multiple columns, separated by comma (`,`).
            These columns will be concatenated to form the cell groups. For the overlapping cells, they will be
            counted multiple times for different groups. So make sure the cell group names in different columns
            are unique.
        cells_order (list): The order of the cells (rows) to show on the plot
        cells_orderby: An expression passed to `dplyr::arrange()` to order the cells (rows) of the plot.
            Only works when `cells-order` is not specified.
            The data frame passed to `dplyr::arrange()` is grouped by `cells_by` before ordering.
            You can have multiple expressions separated by semicolon (`;`). The expessions will be parsed by `rlang::parse_exprs()`.
            4 extra columns were added to the metadata for ordering the rows in the plot:
            * `CloneSize`: The size (number of cells) of clones (identified by `cells_by`)
            * `CloneGroupSize`: The clone size in each group (identified by `group_by`)
            * `CloneClusterSize`: The clone size in each cluster (identified by `seurat_clusters`)
            * `CloneGroupClusterSize`: The clone size in each group and cluster (identified by `group_by` and `seurat_clusters`)
        cells_n (type=int): The max number of groups to show for each cell group identity (row).
            Ignored if `cells_order` is specified.
        subset: An expression to subset the cells, will be passed to `dplyr::filter()` on metadata.
            This will be applied prior to `each`.
        descr: The description of the case, will be shown in the report.
        hm_devpars (ns): The device parameters for the heatmaps.
            - res (type=int): The resolution of the heatmaps.
            - height (type=int): The height of the heatmaps.
            - width (type=int): The width of the heatmaps.
        devpars (ns): The device parameters for the plots of pie charts.
            - res (type=int): The resolution of the plots
            - height (type=int): The height of the plots
            - width (type=int): The width of the plots
        each: The column name in metadata to separate the cells into different plots.
        prefix_each (flag): Whether to prefix the `each` column name to the
            value as the case/section name.
        section: The section to show in the report. This allows different cases to be put in the same section in report.
            Only works when `each` is not specified.
            %(envs_section_each)s
        overlap (list): Plot the overlap of cell groups (values of `cells_by`) in different cases
            under the same section.
            The section must have at least 2 cases, each case should have a single `cells_by` column.
        cases (type=json;order=99): If you have multiple cases, you can specify them here.
            Keys are the names of the cases and values are the options above except `mutaters`.
            If some options are not specified, the options in `envs` will be used.
            If no cases are specified, a default case will be used with case name `DEFAULT`.

    Requires:
        r-seurat:
            - check: {{proc.lang}} -e "library(Seurat)"
        r-dplyr:
            - check: {{proc.lang}} -e "library(dplyr)"
        r-tidyr:
            - check: {{proc.lang}} -e "library(tidyr)"
    """  # noqa: E501
    input = "srtobj:file"
    output = "outdir:dir:{{in.srtobj | stem}}.cells_distribution"
    lang = config.lang.rscript
    envs = {
        "mutaters": {},
        "cluster_orderby": None,
        "group_by": None,
        "group_order": [],
        "cells_by": None,
        "cells_order": [],
        "cells_orderby": None,
        "cells_n": 10,
        "subset": None,
        "descr": None,
        "devpars": {},
        "hm_devpars": {},
        "each": None,
        "prefix_each": True,
        "section": "DEFAULT",
        "overlap": [],
        "cases": {},
    }
    script = "file://../scripts/scrna/CellsDistribution.R"
    plugin_opts = {
        "report": "file://../reports/scrna/CellsDistribution.svelte",
        "report_paging": 8,
    }


@format_placeholder(mutate_helpers_clonesize=MUTATE_HELPERS_CLONESIZE_INDENTED)
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
        mutaters (type=json): The mutaters to mutate the metadata.
            The key-value pairs will be passed the `dplyr::mutate()` to mutate the metadata.
            %(mutate_helpers_clonesize)s

    Requires:
        r-seurat:
            - check: {{proc.lang}} <(echo "library(Seurat)")
        r-tibble:
            - check: {{proc.lang}} <(echo "library(tibble)")
        r-dplyr:
            - check: {{proc.lang}} <(echo "library(dplyr)")
    """  # noqa: E501
    input = "srtobj:file, metafile:file"
    output = "rdsfile:file:{{in.srtobj | stem}}.RDS"
    lang = config.lang.rscript
    envs = {"mutaters": {}}
    script = "file://../scripts/scrna/SeuratMetadataMutater.R"


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


@format_placeholder(
    mutate_helpers_clonesize=MUTATE_HELPERS_CLONESIZE_INDENTED,
    envs_section_each=ENVS_SECTION_EACH_INDENTED,
)
class MarkersFinder(Proc):
    """Find markers between different groups of cells

    When only `group-by` is specified as `"seurat_clusters"` in
    `envs.cases`, the markers will be found for all the clusters.

    You can also find the differentially expressed genes between
    any two groups of cells by setting `group-by` to a different
    column name in metadata. Follow `envs.cases` for more details.

    Input:
        srtobj: The seurat object loaded by `SeuratPreparing`
            If you have your `Seurat` object prepared by yourself, you can also
            use it here, but you should make sure that the object has been processed
            by `PrepSCTFindMarkers` if data is not normalized using `SCTransform`.

    Output:
        outdir: The output directory for the markers

    Envs:
        ncores (type=int): Number of cores to use for parallel computing for some `Seurat` procedures.
            * Used in `future::plan(strategy = "multicore", workers = <ncores>)` to parallelize some Seurat procedures.
            * See also: <https://satijalab.org/seurat/articles/future_vignette.html>
        mutaters (type=json): The mutaters to mutate the metadata
            %(mutate_helpers_clonesize)s
        ident-1: The first group of cells to compare
        ident-2: The second group of cells to compare
            If not provided, the rest of the cells are used for `ident-2`.
        group-by: The column name in metadata to group the cells.
            If only `group-by` is specified, and `ident-1` and `ident-2` are
            not specified, markers will be found for all groups in this column
            in the manner of "group vs rest" comparison.
            `NA` group will be ignored.
        each: The column name in metadata to separate the cells into different
            cases.
        prefix_each (flag): Whether to prefix the `each` column name to the
            value as the case/section name.
        prefix_group (flag): When neither `ident-1` nor `ident-2` is specified,
            should we prefix the group name to the section name?
        dbs (list): The dbs to do enrichment analysis for significant
            markers See below for all libraries.
            <https://maayanlab.cloud/Enrichr/#libraries>
        sigmarkers: An expression passed to `dplyr::filter()` to filter the
            significant markers for enrichment analysis.
            Available variables are `p_val`, `avg_log2FC`, `pct.1`, `pct.2` and
            `p_val_adj`. For example, `"p_val_adj < 0.05 & abs(avg_log2FC) > 1"`
            to select markers with adjusted p-value < 0.05 and absolute log2
            fold change > 1.
        assay: The assay to use.
        volcano_genes (type=auto): The genes to label in the volcano plot if they are
            significant markers.
            If `True`, all significant markers will be labeled. If `False`, no
            genes will be labeled. Otherwise, specify the genes to label.
            It could be either a string with comma separated genes, or a list
            of genes.
        section: The section name for the report. It must not contain colon (`:`).
            Ignored when `each` is not specified and `ident-1` is specified.
            When neither `each` nor `ident-1` is specified, case name will be used
            as section name.
            If `each` is specified, the section name will be constructed from
            `each` and case name.
            %(envs_section_each)s
        subset: An expression to subset the cells for each case.
        rest (ns): Rest arguments for `Seurat::FindMarkers()`.
            Use `-` to replace `.` in the argument name. For example,
            use `min-pct` instead of `min.pct`.
            This only works when `use_presto` is `False`.
            - <more>: See <https://satijalab.org/seurat/reference/findmarkers>
        dotplot (ns): Arguments for `Seurat::DotPlot()`.
            Use `-` to replace `.` in the argument name. For example,
            use `group-bar` instead of `group.bar`.
            Note that `object`, `features`, and `group-by` are already specified
            by this process. So you don't need to specify them here.
            - maxgenes (type=int): The maximum number of genes to plot.
            - devpars (ns): The device parameters for the plots.
                - res (type=int): The resolution of the plots.
                - height (type=int): The height of the plots.
                - width (type=int): The width of the plots.
            - <more>: See <https://satijalab.org/seurat/reference/doheatmap>
        cases (type=json): If you have multiple cases, you can specify them
            here. The keys are the names of the cases and the values are the
            above options except `ncores` and `mutaters`. If some options are
            not specified, the default values specified above will be used.
            If no cases are specified, the default case will be added with
            the default values under `envs` with the name `DEFAULT`.
        overlap_defaults (ns): The default options for overlapping analysis.
            - venn (ns): The options for the Venn diagram.
                Venn diagram can only be plotted for sections with no more than 4 cases.
                - devpars (ns): The device parameters for the plots.
                    - res (type=int): The resolution of the plots.
                    - height (type=int): The height of the plots.
                    - width (type=int): The width of the plots.
            - upset (ns): The options for the UpSet plot.
                - devpars (ns): The device parameters for the plots.
                    - res (type=int): The resolution of the plots.
                    - height (type=int): The height of the plots.
                    - width (type=int): The width of the plots.
        overlap (json): The sections to do overlaping analysis, including
            Venn diagram and UpSet plot. The Venn diagram and UpSet plot
            will be plotted for the overlapping of significant markers between
            different cases.
            The keys of this option are the names of the sections. The values are
            a dict of options with keys `venn` and `upset`, values will
            be inherited from `envs.overlap_defaults`, recursively.
            You can set `envs.overlap.<section>.venn` to `False`/`None` to disable
            the Venn diagram for the section.
            It works when `each` is specified. In such a case, the sections will be
            the case names.
            This does not work for the cases where `ident-1` is not specified. In case
            you want to do such analysis for those cases, you should enumerate the
            idents in different cases and specify them here.
        cache (type=auto): Where to cache to `FindAllMarkers` results.
            If `True`, cache to `outdir` of the job. If `False`, don't cache.
            Otherwise, specify the directory to cache to.
            Only works when `use_presto` is `False` (presto works fast enough).
    """  # noqa: E501
    input = "srtobj:file"
    output = "outdir:dir:{{in.srtobj | stem0}}.markers"
    lang = config.lang.rscript
    envs = {
        "ncores": config.misc.ncores,
        "mutaters": {},
        "ident-1": None,
        "ident-2": None,
        "group-by": "seurat_clusters",
        "each": None,
        "prefix_each": True,
        "prefix_group": True,
        "section": "DEFAULT",
        "assay": None,
        "subset": None,
        "rest": {},
        "dbs": ["KEGG_2021_Human", "MSigDB_Hallmark_2020"],
        "sigmarkers": "p_val_adj < 0.05",
        "volcano_genes": True,
        "dotplot": {"maxgenes": 20},
        "cases": {},
        "overlap_defaults": {
            "venn": {"devpars": {"res": 100, "height": 600, "width": 1000}},
            "upset": {"devpars": {"res": 100, "height": 600, "width": 800}},
        },
        "overlap": {},
        "cache": config.path.tmpdir,
    }
    order = 5
    script = "file://../scripts/scrna/MarkersFinder.R"
    plugin_opts = {
        "report": "file://../reports/scrna/MarkersFinder.svelte",
        "report_paging": 8,
    }


class TopExpressingGenes(Proc):
    """Find the top expressing genes in each cluster

    Input:
        srtobj: The seurat object in RDS format

    Output:
        outdir: The output directory for the tables and plots

    Envs:
        mutaters (type=json): The mutaters to mutate the metadata
        ident: The group of cells to find the top expressing genes.
            The cells will be selected by the `group-by` column with this
            `ident` value in metadata.
            If not provided, the top expressing genes will be found for all
            groups of cells in the `group-by` column.
        group-by: The column name in metadata to group the cells.
        each: The column name in metadata to separate the cells into different
            cases.
            When specified, `ident` must be specified
        prefix_each (flag): Whether to prefix the `each` column name to the
            value as the case/section name.
        section: The section name for the report.
            Worked only when `each` is not specified and `ident` is specified.
            Otherwise, the section name will be constructed from `each` and
            `group-by`.
            If `DEFAULT`, and it's the only section, it not included in the
            case/section names.
        dbs (list): The dbs to do enrichment analysis for significant
            markers See below for all libraries.
            <https://maayanlab.cloud/Enrichr/#libraries>
        n (type=int): The number of top expressing genes to find.
        subset: An expression to subset the cells for each case.
        cases (type=json): If you have multiple cases, you can specify them
            here. The keys are the names of the cases and the values are the
            above options except `mutaters`. If some options are
            not specified, the default values specified above will be used.
            If no cases are specified, the default case will be added with
            the default values under `envs` with the name `DEFAULT`.
    """
    input = "srtobj:file"
    output = "outdir:dir:{{in.srtobj | stem}}.top_expressing_genes"
    lang = config.lang.rscript
    script = "file://../scripts/scrna/TopExpressingGenes.R"
    envs = {
        "mutaters": {},
        "ident": None,
        "group-by": "seurat_clusters",
        "each": None,
        "prefix_each": True,
        "section": "DEFAULT",
        "dbs": ["KEGG_2021_Human", "MSigDB_Hallmark_2020"],
        "n": 250,
        "subset": None,
        "cases": {},
    }
    plugin_opts = {
        "report": "file://../reports/scrna/TopExpressingGenes.svelte",
        "report_paging": 8,
    }


class ExprImputation(Proc):
    """This process imputes the dropout values in scRNA-seq data.

    It takes the Seurat object as input and outputs the Seurat object with
    imputed expression data.

    Reference:
    - [Linderman, George C., Jun Zhao, and Yuval Kluger. "Zero-preserving imputation of scRNA-seq data using low-rank approximation." BioRxiv (2018): 397588.](https://www.nature.com/articles/s41467-021-27729-z)
    - [Li, Wei Vivian, and Jingyi Jessica Li. "An accurate and robust imputation method scImpute for single-cell RNA-seq data." Nature communications 9.1 (2018): 997.](https://www.nature.com/articles/s41467-018-03405-7)
    - [Dijk, David van, et al. "MAGIC: A diffusion-based imputation method reveals gene-gene interactions in single-cell RNA-sequencing data." BioRxiv (2017): 111591.](https://www.cell.com/cell/abstract/S0092-8674(18)30724-4)

    Input:
        infile: The input file in RDS format of Seurat object

    Output:
        outfile: The output file in RDS format of Seurat object
            Note that with rmagic and alra, the original default assay will be
            renamed to `RAW` and the imputed RNA assay will be
            renamed to `RNA` and set as default assay.

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
        alra_args (type=json): The arguments for `RunALRA()`

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
    """  # noqa: E501
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
    script = "file://../scripts/scrna/ExprImputation.R"


class SCImpute(Proc):
    """Impute the dropout values in scRNA-seq data.

    Deprecated. Use `ExprImputation` instead.

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


class SeuratTo10X(Proc):
    """Write a Seurat object to 10X format

    using `write10xCounts` from `DropletUtils`

    Input:
        srtobj: The seurat object in RDS

    Output:
        outdir: The output directory.
            When `envs.split_by` is specified, the subdirectories will be
            created for each distinct value of the column.
            Otherwise, the matrices will be written to the output directory.

    Envs:
        version: The version of 10X format
    """
    input = "srtobj:file"
    output = "outdir:dir:{{in.srtobj | stem}}"
    envs = {"version": "3", "split_by": None}
    lang = config.lang.rscript
    script = "file://../scripts/scrna/SeuratTo10X.R"


@format_placeholder(
    mutate_helpers_clonesize=MUTATE_HELPERS_CLONESIZE_INDENTED,
    envs_section_each=ENVS_SECTION_EACH_INDENTED,
)
class ScFGSEA(Proc):
    """Gene set enrichment analysis for cells in different groups using `fgsea`

    This process allows us to do Gene Set Enrichment Analysis (GSEA) on the expression data,
    but based on variaties of grouping, including the from the meta data and the
    scTCR-seq data as well.

    The GSEA is done using the
    [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html) package,
    which allows to quickly and accurately calculate arbitrarily low GSEA P-values
    for a collection of gene sets.
    The fgsea package is based on the fast algorithm for preranked GSEA described in
    [Subramanian et al. 2005](https://www.pnas.org/content/102/43/15545).

    For each case, the process will generate a table with the enrichment scores for
    each gene set, and GSEA plots for the top gene sets.

    Input:
        srtobj: The seurat object in RDS format

    Output:
        outdir: The output directory for the results

    Envs:
        ncores (type=int): Number of cores for parallelization
            Passed to `nproc` of `fgseaMultilevel()`.
        mutaters (type=json): The mutaters to mutate the metadata.
            The key-value pairs will be passed the `dplyr::mutate()` to mutate the metadata.
            %(mutate_helpers_clonesize)s

        group-by: The column name in metadata to group the cells.
        ident-1: The first group of cells to compare
        ident-2: The second group of cells to compare, if not provided, the rest of the cells that are not `NA`s in `group-by` column are used for `ident-2`.
        each: The column name in metadata to separate the cells into different subsets to do the analysis.
        prefix_each (flag): Whether to prefix the `each` column name to the values as the case/section name.
        subset: An expression to subset the cells.
        section: The section name for the report. Worked only when `each` is not specified. Otherwise, the section name will be constructed from `each` and its value.
            This allows different cases to be put into the same section in the report.
            %(envs_section_each)s
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
            If no cases are specified, the default case will be added with the name `DEFAULT`.

    Requires:
        bioconductor-fgsea:
            - check: {{proc.lang}} -e "library(fgsea)"
        r-seurat:
            - check: {{proc.lang}} -e "library(seurat)"
    """  # noqa: E501
    input = "srtobj:file"
    output = "outdir:dir:{{(in.casefile or in.srtobj) | stem0}}.fgsea"
    lang = config.lang.rscript
    envs = {
        "mutaters": {},
        "ncores": config.misc.ncores,
        "group-by": None,
        "ident-1": None,
        "ident-2": None,
        "each": None,
        "prefix_each": True,
        "subset": None,
        "section": "DEFAULT",
        "gmtfile": "",
        "method": "s2n",
        "top": 20,
        "minsize": 10,
        "maxsize": 100,
        "eps": 0,
        "rest": {},
        "cases": {},
    }
    script = "file://../scripts/scrna/ScFGSEA.R"
    plugin_opts = {
        "report": "file://../reports/scrna/ScFGSEA.svelte",
        "report_paging": 8,
    }


class CellTypeAnnotation(Proc):
    """Annotate the cell clusters. Currently, four ways are supported:

    1. Pass the cell type annotation directly
    2. Use [`ScType`](https://github.com/IanevskiAleksandr/sc-type)
    3. Use [`scCATCH`](https://github.com/ZJUFanLab/scCATCH)
    4. Use [`hitype`](https://github.com/pwwang/hitype)

    The annotated cell types will replace the original `seurat_clusters` column in the metadata,
    so that the downstream processes will use the annotated cell types.

    The old `seurat_clusters` column will be renamed to `seurat_clusters_id`.

    If you are using `ScType`, `scCATCH`, or `hitype`, a text file containing the mapping from
    the old `seurat_clusters` to the new cell types will be generated and saved to
    `cluster2celltype.tsv` under `<workdir>/<pipline_name>/CellTypeAnnotation/0/output/`.

    Examples:
        ```toml
        [CellTypeAnnotation.envs]
        tool = "direct"
        cell_types = ["CellType1", "CellType2", "-", "CellType4"]
        ```

        The cell types will be assigned as:

        ```
        0 -> CellType1
        1 -> CellType2
        2 -> 2
        3 -> CellType4
        ```

    Input:
        sobjfile: The seurat object

    Output:
        outfile: The rds file of seurat object with cell type annotated

    Envs:
        tool (choice): The tool to use for cell type annotation.
            - sctype: Use `scType` to annotate cell types.
                See <https://github.com/IanevskiAleksandr/sc-type>
            - hitype: Use `hitype` to annotate cell types.
                See <https://github.com/pwwang/hitype>
            - sccatch: Use `scCATCH` to annotate cell types.
                See <https://github.com/ZJUFanLab/scCATCH>
            - celltypist: Use `celltypist` to annotate cell types.
                See <https://github.com/Teichlab/celltypist>
            - direct: Directly assign cell types
        sctype_tissue: The tissue to use for `sctype`.
            Avaiable tissues should be the first column (`tissueType`) of `sctype_db`.
            If not specified, all rows in `sctype_db` will be used.
        sctype_db: The database to use for sctype.
            Check examples at <https://github.com/IanevskiAleksandr/sc-type/blob/master/ScTypeDB_full.xlsx>
        hitype_tissue: The tissue to use for `hitype`.
            Avaiable tissues should be the first column (`tissueType`) of `hitype_db`.
            If not specified, all rows in `hitype_db` will be used.
        hitype_db: The database to use for hitype.
            Compatible with `sctype_db`.
            See also <https://pwwang.github.io/hitype/articles/prepare-gene-sets.html>
            You can also use built-in databases, including `hitypedb_short`, `hitypedb_full`, and `hitypedb_pbmc3k`.
        cell_types (list): The cell types to use for direct annotation.
            You can use `"-"` or `""` as the placeholder for the clusters that
            you want to keep the original cell types (`seurat_clusters`).
            If the length of `cell_types` is shorter than the number of
            clusters, the remaining clusters will be kept as the original cell
            types.
            You can also use `NA` to remove the clusters from downstream analysis. This
            only works when `envs.newcol` is not specified.

            /// Note
            If `tool` is `direct` and `cell_types` is not specified or an empty list,
            the original cell types will be kept and nothing will be changed.
            ///

        sccatch_args (ns): The arguments for `scCATCH::findmarkergene()` if `tool` is `sccatch`.
            - species: The specie of cells.
            - cancer: If the sample is from cancer tissue, then the cancer type may be defined.
            - tissue: Tissue origin of cells must be defined.
            - marker: The marker genes for cell type identification.
            - if_use_custom_marker (flag): Whether to use custom marker genes. If `True`, no `species`, `cancer`, and `tissue` are needed.
            - <more>: Other arguments for [`scCATCH::findmarkergene()`](https://rdrr.io/cran/scCATCH/man/findmarkergene.html).
                You can pass an RDS file to `sccatch_args.marker` to work as custom marker. If so,
                `if_use_custom_marker` will be set to `TRUE` automatically.
        celltypist_args (ns): The arguments for `celltypist::celltypist()` if `tool` is `celltypist`.
            - model: The path to model file.
            - python: The python path where celltypist is installed.
            - majority_voting: When true, it refines cell identities within local subclusters after an over-clustering approach
                at the cost of increased runtime.
            - over_clustering (type=auto): The column name in metadata to use as clusters for majority voting.
                Set to `False` to disable over-clustering.
            - assay: When converting a Seurat object to AnnData, the assay to use.
                If input is h5seurat, this defaults to RNA.
                If input is Seurat object in RDS, this defaults to the default assay.
        newcol: The new column name to store the cell types.
            If not specified, the `seurat_clusters` column will be overwritten.
            If specified, the original `seurat_clusters` column will be kept and `Idents` will be kept as the original `seurat_clusters`.
        outtype (choice): The output file type. Currently only works for `celltypist`.
            An RDS file will be generated for other tools.
            - input: Use the same file type as the input.
            - rds: Use RDS file.
            - h5seurat: Use h5seurat file.
            - h5ad: Use AnnData file.

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
    """  # noqa: E501
    input = "sobjfile:file"
    output = (
        "outfile:file:"
        "{{in.sobjfile | stem}}.annotated."
        "{{- ext0(in.sobjfile) if envs.outtype == 'input' else envs.outtype -}}"
    )
    lang = config.lang.rscript
    envs = {
        "tool": "hitype",
        "sctype_tissue": None,
        "sctype_db": config.ref.sctype_db,
        "cell_types": [],
        "sccatch_args": {
            "species": None,
            "cancer": "Normal",
            "tissue": None,
            "marker": None,
            "if_use_custom_marker": False,
        },
        "hitype_tissue": None,
        "hitype_db": None,
        "celltypist_args": {
            "model": None,
            "python": config.lang.python,
            "majority_voting": True,
            "over_clustering": "seurat_clusters",
            "assay": None,
        },
        "newcol": None,
        "outtype": "input",
    }
    script = "file://../scripts/scrna/CellTypeAnnotation.R"


class SeuratMap2Ref(Proc):
    """Map the seurat object to reference

    See: <https://satijalab.org/seurat/articles/integration_mapping.html>
    and <https://satijalab.org/seurat/articles/multimodal_reference_mapping.html>

    Input:
        sobjfile: The seurat object

    Output:
        outfile: The rds file of seurat object with cell type annotated.
            Note that the reduction name will be `ref.umap` for the mapping.
            To visualize the mapping, you should use `ref.umap` as the reduction name.

    Envs:
        ncores (type=int;order=-100): Number of cores to use.
            When `split_by` is used, this will be the number of cores for each object to map to the reference.
            When `split_by` is not used, this is used in `future::plan(strategy = "multicore", workers = <ncores>)`
            to parallelize some Seurat procedures.
            See also: <https://satijalab.org/seurat/archive/v3.0/future_vignette.html>
        mutaters (type=json): The mutaters to mutate the metadata.
            This is helpful when we want to create new columns for `split_by`.
        use: A column name of metadata from the reference
            (e.g. `celltype.l1`, `celltype.l2`) to transfer to the query as the
            cell types (ident) for downstream analysis. This field is required.
            If you want to transfer multiple columns, you can use
            `envs.MapQuery.refdata`.
        ident: The name of the ident for query transferred from `envs.use` of the reference.
        ref: The reference seurat object file.
            Either an RDS file or a h5seurat file that can be loaded by
            `Seurat::LoadH5Seurat()`.
            The file type is determined by the extension. `.rds` or `.RDS` for
            RDS file, `.h5seurat` or `.h5` for h5seurat file.
        refnorm (choice): Normalization method the reference used. The same method will be used for the query.
            - NormalizeData: Using [`NormalizeData`](https://satijalab.org/seurat/reference/normalizedata).
            - SCTransform: Using [`SCTransform`](https://satijalab.org/seurat/reference/sctransform).
            - auto: Automatically detect the normalization method.
                If the default assay of reference is `SCT`, then `SCTransform` will be used.
        split_by: The column name in metadata to split the query into multiple objects.
            This helps when the original query is too large to process.
        SCTransform (ns): Arguments for [`SCTransform()`](https://satijalab.org/seurat/reference/sctransform)
            - do-correct-umi (flag): Place corrected UMI matrix in assay counts layer?
            - do-scale (flag): Whether to scale residuals to have unit variance?
            - do-center (flag): Whether to center residuals to have mean zero?
            - <more>: See <https://satijalab.org/seurat/reference/sctransform>.
                Note that the hyphen (`-`) will be transformed into `.` for the keys.
        NormalizeData (ns): Arguments for [`NormalizeData()`](https://satijalab.org/seurat/reference/normalizedata)
            - normalization-method: Normalization method.
            - <more>: See <https://satijalab.org/seurat/reference/normalizedata>.
                Note that the hyphen (`-`) will be transformed into `.` for the keys.
        FindTransferAnchors (ns): Arguments for [`FindTransferAnchors()`](https://satijalab.org/seurat/reference/findtransferanchors)
            - normalization-method (choice): Name of normalization method used.
                - LogNormalize: Log-normalize the data matrix
                - SCT: Scale data using the SCTransform method
                - auto: Automatically detect the normalization method.
                    See `envs.refnorm`.
            - reference-reduction: Name of dimensional reduction to use from the reference if running the pcaproject workflow.
                Optionally enables reuse of precomputed reference dimensional reduction.
            - <more>: See <https://satijalab.org/seurat/reference/findtransferanchors>.
                Note that the hyphen (`-`) will be transformed into `.` for the keys.
        MapQuery (ns): Arguments for [`MapQuery()`](https://satijalab.org/seurat/reference/mapquery)
            - reference-reduction: Name of reduction to use from the reference for neighbor finding
            - reduction-model: `DimReduc` object that contains the umap model.
            - refdata (type=json): Extra data to transfer from the reference to the query.
            - <more>: See <https://satijalab.org/seurat/reference/mapquery>.
                Note that the hyphen (`-`) will be transformed into `.` for the keys.
        MappingScore (ns): Arguments for [`MappingScore()`](https://satijalab.org/seurat/reference/mappingscore)
            - <more>: See <https://satijalab.org/seurat/reference/mappingscore>.
                Note that the hyphen (`-`) will be transformed into `.` for the keys.

    Requires:
        r-seurat:
            - check: {{proc.lang}} -e "library(Seurat)"
    """  # noqa: E501
    input = "sobjfile:file"
    output = "outfile:file:{{in.sobjfile | stem}}.RDS"
    lang = config.lang.rscript
    envs = {
        "ncores": config.misc.ncores,
        "use": None,
        "ident": "seurat_clusters",
        "mutaters": {},
        "ref": None,
        "refnorm": "auto",
        "split_by": None,
        "SCTransform": {
            "do-correct-umi": False,
            "do-scale": False,
            "do-center": True,
        },
        "NormalizeData": {
            "normalization-method": "LogNormalize",
        },
        "FindTransferAnchors": {
            "reference-reduction": "spca",
        },
        "MapQuery": {
            "reference-reduction": "spca",
            "reduction-model": "wnn.umap",
            "refdata": {
                # "celltype-l1": "celltype.l1",
                # "celltype-l2": "celltype.l2",
                # "predicted_ADT": "ADT",
            }
        },
        "MappingScore": {"ndim": 30},
    }
    script = "file://../scripts/scrna/SeuratMap2Ref.R"
    plugin_opts = {"report": "file://../reports/scrna/SeuratMap2Ref.svelte"}


class RadarPlots(Proc):
    """Radar plots for cell proportion in different clusters.

    This process generates the radar plots for the clusters of T cells.
    It explores the proportion of cells in different groups (e.g. Tumor vs Blood)
    in different T-cell clusters.

    Examples:
        Let's say we have a metadata like this:

        | Cell | Source | Timepoint | seurat_clusters |
        | ---- | ------ | --------- | --------------- |
        | A    | Blood  | Pre       | 0               |
        | B    | Blood  | Pre       | 0               |
        | C    | Blood  | Post      | 1               |
        | D    | Blood  | Post      | 1               |
        | E    | Tumor  | Pre       | 2               |
        | F    | Tumor  | Pre       | 2               |
        | G    | Tumor  | Post      | 3               |
        | H    | Tumor  | Post      | 3               |

        With configurations:

        ```toml
        [RadarPlots.envs]
        by = "Source"
        ```

        Then we will have a radar plots like this:

        ![Radar plots](https://pwwang.github.io/immunopipe/latest/processes/images/RadarPlots-default.png)

        We can use `each` to separate the cells into different cases:

        ```toml
        [RadarPlots.envs]
        by = "Source"
        each = "Timepoint"
        ```

        Then we will have two radar plots, one for `Pre` and one for `Post`:

        ![Radar plots](https://pwwang.github.io/immunopipe/latest/processes/images/RadarPlots-each.png)

        Using `cluster_order` to change the order of the clusters and show only the first 3 clusters:

        ```toml
        [RadarPlots.envs]
        by = "Source"
        cluster_order = ["2", "0", "1"]
        breaks = [0, 50, 100]  # also change the breaks
        ```

        ![Radar plots cluster_order](https://pwwang.github.io/immunopipe/latest/processes/images/RadarPlots-cluster_order.png)


        /// Attention
        All the plots used in the examples are just for demonstration purpose. The real plots will have different appearance.
        ///

    Input:
        srtobj: The seurat object in RDS format

    Output:
        outdir: The output directory for the plots

    Envs:
        mutaters (type=json): Mutaters to mutate the metadata of the
            seurat object. Keys are the column names and values are the
            expressions to mutate the columns. These new columns will be
            used to define your cases.
        by: Which column to use to separate the cells in different groups.
            `NA`s will be ignored. For example, If you have a column named `Source`
            that marks the source of the cells, and you want to separate the cells
            into `Tumor` and `Blood` groups, you can set `by` to `Source`.
            The there will be two curves in the radar plot, one for `Tumor` and
            one for `Blood`.
        each: A column with values to separate all cells in different cases
            When specified, the case will be expanded to multiple cases for
            each value in the column.
            If specified, `section` will be ignored, and the case name will
            be used as the section name.
        prefix_each (flag): Whether to prefix the `each` column name to the values as the
            case/section name.
        breakdown: An additional column with groups to break down the cells
            distribution in each cluster. For example, if you want to see the
            distribution of the cells in each cluster in different samples. In
            this case, you should have multiple values in each `by`. These values
            won't be plotted in the radar plot, but a barplot will be generated
            with the mean value of each group and the error bar.
        test (choice): The test to use to calculate the p values.
            If there are more than 2 groups in `by`, the p values will be calculated
            pairwise group by group. Only works when `breakdown` is specified and
            `by` has 2 groups or more.
            - wilcox: Wilcoxon rank sum test
            - t: T test
            - none: No test will be performed
        order (list): The order of the values in `by`. You can also limit
            (filter) the values we have in `by`. For example, if column `Source`
            has values `Tumor`, `Blood`, `Spleen`, and you only want to plot
            `Tumor` and `Blood`, you can set `order` to `["Tumor", "Blood"]`.
            This will also have `Tumor` as the first item in the legend and `Blood`
            as the second item.
        colors: The colors for the groups in `by`. If not specified,
            the default colors will be used.
            Multiple colors can be separated by comma (`,`).
            You can specify `biopipen` to use the `biopipen` palette.
        ident: The column name of the cluster information.
        cluster_order (list): The order of the clusters.
            You may also use it to filter the clusters. If not given,
            all clusters will be used.
            If the cluster names are integers, use them directly for the order,
            even though a prefix `Cluster` is added on the plot.
        breaks (list;itype=int): breaks of the radar plots, from 0 to 100.
            If not given, the breaks will be calculated automatically.
        direction (choice): Direction to calculate the percentages.
            - inter-cluster: the percentage of the cells in all groups
                in each cluster (percentage adds up to 1 for each cluster).
            - intra-cluster: the percentage of the cells in all clusters.
                (percentage adds up to 1 for each group).
        section: If you want to put multiple cases into a same section
            in the report, you can set this option to the name of the section.
            Only used in the report.
        subset: The subset of the cells to do the analysis.
        bar_devpars (ns): The parameters for `png()` for the barplot
            - res (type=int): The resolution of the plot
            - height (type=int): The height of the plot
            - width (type=int): The width of the plot
        devpars (ns): The parameters for `png()`
            - res (type=int): The resolution of the plot
            - height (type=int): The height of the plot
            - width (type=int): The width of the plot
        cases (type=json): The cases for the multiple radar plots.
            Keys are the names of the cases and values are the arguments for
            the plots (`each`, `by`, `order`, `breaks`, `direction`,
            `ident`, `cluster_order` and `devpars`).
            If not cases are given, a default case will be used, with the
            key `DEFAULT`.
            The keys must be valid string as part of the file name.
    """  # noqa: E501
    input = "srtobj:file"
    output = "outdir:dir:{{in.srtobj | stem}}.radar_plots"
    lang = config.lang.rscript
    script = "file://../scripts/scrna/RadarPlots.R"
    envs = {
        "mutaters": {},
        "by": None,
        "each": None,
        "prefix_each": True,
        "order": None,
        "colors": None,
        "ident": "seurat_clusters",
        "cluster_order": [],
        "breakdown": None,
        "test": "wilcox",
        "breaks": [],
        "direction": "intra-cluster",
        "section": "DEFAULT",
        "subset": None,
        "bar_devpars": {
            "res": 100,
            "width": 1200,
            "height": 800,
        },
        "devpars": {
            "res": 100,
            "width": 1200,
            "height": 1000,
        },
        "cases": {},
    }
    plugin_opts = {
        "report": "file://../reports/scrna/RadarPlots.svelte",
    }


@format_placeholder(
    mutate_helpers_clonesize=MUTATE_HELPERS_CLONESIZE_INDENTED,
    envs_section_each=ENVS_SECTION_EACH_INDENTED,
)
class MetaMarkers(Proc):
    """Find markers between three or more groups of cells, using one-way ANOVA
    or Kruskal-Wallis test.

    Sometimes, you may want to find the markers for cells from more than 2 groups.
    In this case, you can use this process to find the markers for the groups and
    do enrichment analysis for the markers. Each marker is examined using either
    one-way ANOVA or Kruskal-Wallis test.
    The p values are adjusted using the specified method. The significant markers
    are then used for enrichment analysis using
    [enrichr](https://maayanlab.cloud/Enrichr/) api.

    Other than the markers and the enrichment analysis as outputs, this process also
    generates violin plots for the top 10 markers.

    Input:
        srtobj: The seurat object loaded by `SeuratPreparing`

    Output:
        outdir: The output directory for the markers

    Envs:
        ncores (type=int): Number of cores to use to parallelize for genes
        mutaters (type=json): The mutaters to mutate the metadata
            The key-value pairs will be passed the `dplyr::mutate()` to mutate the metadata.
            %(mutate_helpers_clonesize)s

        group-by: The column name in metadata to group the cells.
            If only `group-by` is specified, and `idents` are
            not specified, markers will be found for all groups in this column.
            `NA` group will be ignored.
        idents: The groups of cells to compare, values should be in the `group-by` column.
        each: The column name in metadata to separate the cells into different cases.
        prefix_each (flag): Whether to add the `each` value as prefix to the case name.
        dbs (list): The dbs to do enrichment analysis for significant
            markers See below for all libraries.
            <https://maayanlab.cloud/Enrichr/#libraries>
        subset: The subset of the cells to do the analysis.
            An expression passed to `dplyr::filter()`.
        p_adjust (choice): The method to adjust the p values, which can be used to filter the significant markers.
            See also <https://rdrr.io/r/stats/p.adjust.html>
            - holm: Holm-Bonferroni method
            - hochberg: Hochberg method
            - hommel: Hommel method
            - bonferroni: Bonferroni method
            - BH: Benjamini-Hochberg method
            - BY: Benjamini-Yekutieli method
            - fdr: FDR method of Benjamini-Hochberg
            - none: No adjustment
        sigmarkers: An expression passed to `dplyr::filter()` to filter the
            significant markers for enrichment analysis. The default is `p.value < 0.05`.
            If `method = 'anova'`, the variables that can be used for filtering are:
            `sumsq`, `meansq`, `statistic`, `p.value` and `p_adjust`.
            If `method = 'kruskal'`, the variables that can be used for filtering are:
            `statistic`, `p.value` and `p_adjust`.
        section: The section name for the report.
            Worked only when `each` is not specified.
            Otherwise, the section name will be constructed from `each` and `group-by`.
            If `DEFAULT`, and it's the only section, it not included in the case/section names.
            %(envs_section_each)s
        method (choice): The method for the test.
            - anova: One-way ANOVA
            - kruskal: Kruskal-Wallis test
        cases (type=json): If you have multiple cases, you can specify them
            here. The keys are the names of the cases and the values are the
            above options except `ncores` and `mutaters`. If some options are
            not specified, the default values specified above will be used.
            If no cases are specified, the default case will be added with
            the default values under `envs` with the name `DEFAULT`.
    """  # noqa: E501
    input = "srtobj:file"
    output = "outdir:dir:{{in.srtobj | stem}}.meta_markers"
    lang = config.lang.rscript
    script = "file://../scripts/scrna/MetaMarkers.R"
    envs = {
        "ncores": config.misc.ncores,
        "mutaters": {},
        "group-by": None,
        "idents": None,
        "each": None,
        "subset": None,
        "prefix_each": True,
        "p_adjust": "BH",
        "dbs": ["KEGG_2021_Human", "MSigDB_Hallmark_2020"],
        "sigmarkers": "p_adjust < 0.05",
        "section": "DEFAULT",
        "method": "anova",
        "cases": {},
    }
    plugin_opts = {
        "report": "file://../reports/scrna/MetaMarkers.svelte",
        "report_paging": 8,
        "poplog_max": 15,
    }


class Seurat2AnnData(Proc):
    """Convert seurat object to AnnData

    Input:
        sobjfile: The seurat object file, in RDS or h5seurat format

    Output:
        outfile: The AnnData file

    Envs:
        assay: The assay to use for AnnData.
            If not specified, the default assay will be used.
    """
    input = "sobjfile:file"
    output = "outfile:file:{{in.sobjfile | stem}}.h5ad"
    lang = config.lang.rscript
    script = "file://../scripts/scrna/Seurat2AnnData.R"
    envs = {"assay": None}


class AnnData2Seurat(Proc):
    """Convert AnnData to seurat object

    Input:
        adfile: The AnnData file

    Output:
        outfile: The seurat object file in RDS format

    Envs:
        assay: The assay to use to convert to seurat object.
        outtype (choice): The output file type.
            - rds: RDS file
            - h5seurat: h5seurat file
        dotplot_check (type=auto): Whether to do a check with `Seurat::DotPlot`
            to see if the conversion is successful.
            Set to `False` to disable the check.
            If `True`, top 10 variable genes will be used for the check.
            You can give a list of genes or a string of genes with comma (`,`) separated
            to use for the check.
            Only works for `outtype = 'rds'`.
    """
    input = "adfile:file"
    output = "outfile:file:{{in.adfile | stem}}.RDS"
    lang = config.lang.rscript
    envs = {"outtype": "rds", "assay": "RNA", "dotplot_check": True}
    script = "file://../scripts/scrna/AnnData2Seurat.R"
