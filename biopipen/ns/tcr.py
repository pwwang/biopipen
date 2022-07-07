"""Tools to analyze single-cell TCR sequencing data"""

from ..core.defaults import SCRIPT_DIR
from ..core.proc import Proc
from ..core.config import config


class ImmunarchLoading(Proc):
    """Immuarch - Loading data

    Build based on immunarch 0.6.7
    See https://immunarch.com/articles/v2_data.html for supported data formats
    Currently only 10x data format is supported

    Library `dplyr` is also required to manipulate the meta data.

    Input:
        metafile: The meta data of the samples
            A tab-delimited file
            Two columns are required:
            - `Sample` to specify the sample names.
            - `TCRDir` to assign the path of the data to the samples,
            and this column will be excluded as metadata.
            Immunarch is able to fetch the sample names from the names of
            the target files. However, 10x data yields result like
            `filtered_contig_annotations.csv`, which doesn't have any name
            information.

    Output:
        rdsfile: The RDS file with the data and metadata
        metatxt: The meta data of the cells, used to attach to the Seurat object

    Envs:
        tmpdir: The temporary directory to link all data files.
        prefix: The prefix to the barcodes. You can use placeholder like
            `{Sample}_` to use the meta data from the immunarch object
        mode: Either "single" for single chain data or "paired" for
            paired chain data. For `single`, only TRB chain will be kept
            at `immdata$data`, information for other chains will be
            saved at `immdata$tra` and `immdata$multi`.
        metacols: The columns to be exported from the metatxt.
    """

    input = "metafile:file"
    output = [
        "rdsfile:file:{{in.metafile | stem}}.immunarch.RDS",
        "metatxt:file:{{in.metafile | stem}}.tcr.txt",
    ]
    lang = config.lang.rscript
    envs = {
        "tmpdir": config.path.tmpdir,
        "prefix": "{Sample}_",
        "mode": "single",
        "metacols": ["Clones", "Proportion", "CDR3.aa"],
    }
    script = "file://../scripts/tcr/ImmunarchLoading.R"


class ImmunarchFilter(Proc):
    """Immunarch - Filter data

    See https://immunarch.com/articles/web_only/repFilter_v3.html

    Input:
        immdata: The data loaded by `immunarch::repLoad()`
        filterfile: A config file in TOML.
            A dict of configurations with keys as the names of the group and
            values dicts with following keys.
            See `envs.filters`


    Output:
        outfile: The filtered `immdata`
        groupfile: Also a group file with rownames as cells and column names as
            each of the keys in `in.filterfile` or `envs.filters`. The values
            will be subkeys of the dicts in `in.filterfile` or `envs.filters`.

    Envs:
        filters: The filters to filter the data
            You can have multiple cases (groups), the names will be the keys of
            this dict, values are also dicts with keys the methods supported by
            `immunarch::repFilter()`.
            There is one more method `by.count` supported to filter the
            count matrix. For `by.meta`, `by.repertoire`, `by.rep`,
            `by.clonotype` or `by.col` the values will be passed to
            `.query` of `repFilter()`.
            You can also use the helper functions provided by `immunarch`,
            including `morethan`, `lessthan`, `include`, `exclude` and
            `interval`. If these functions are not used, `include(value)` will
            be used by default.
            For `by.count`, the value of `filter` will be passed to
            `dplyr::filter()` to filter the count matrix.
            You can also specify `ORDER` to define the filtration order, which
            defaults to 0, higher `ORDER` gets later executed.
            Each subkey/subgroup must be exclusive
            For example:
            >>> {
            >>>   "name": "BM_Post_Clones",
            >>>   "filters" {
            >>>     "Top_20": {
            >>>       "SAVE": True,  # Save the filtered data to immdata
            >>>       "by.meta": {"Source": "BM", "Status": "Post"},
            >>>       "by.count": {
            >>>         "ORDER": 1, "filter": "TOTAL %%in%% TOTAL[1:20]"
            >>>        }
            >>>     },
            >>>     "Rest": {
            >>>       "by.meta": {"Source": "BM", "Status": "Post"},
            >>>       "by.count": {
            >>>         "ORDER": 1, "filter": "!TOTAL %%in%% TOTAL[1:20]"
            >>>        }
            >>>   }
            >>> }

        prefix: The prefix will be added to the cells in the output file
            Placeholders like `{Sample}_` can be used to from the meta data
        metacols: The extra columns to be exported to the group file.
    """
    input = "immdata:file, filterfile:file"
    output = """
        outfile:file:{{in.immdata | stem}}.RDS,
        groupfile:file:{% if in.filterfile -%}
            {{- in.filterfile | toml_load | attr: "name" | append: ".txt" -}}
        {%- else -%}
            {{- envs.filters | attr: "name" | append: ".txt" -}}
        {%- endif -%}
    """
    envs = {
        "prefix": "{Sample}_",
        "filters": {},
        "metacols": ["Clones", "Proportion", "CDR3.aa"],
    }
    lang = config.lang.rscript
    script = "file://../scripts/tcr/ImmunarchFilter.R"


class Immunarch(Proc):
    """Exploration of Single-cell and Bulk T-cell/Antibody Immune Repertoires

    See https://immunarch.com/articles/web_only/v3_basic_analysis.html

    Input:
        immdata: The data loaded by `immunarch::repLoad()`

    Output:
        outdir: The output directory

    Envs:
        volume_by: Groupings to show clonotype volume (sizes)
            Multiple groups supported, for example:
            `volume_by = {0: "Status", 1: ["Status", "Sex"]}`
            Or label the groups:
            `volume_by = {"Status": "Status", "Status_Sex": ["Status", "Sex"]}`
            If a list or a single variable is given, it will be changed
            into `{"Status": "Status"}`
        len_by: Groupings to show CDR3 length of both aa and nt
        count_by: Groupings to show clonotype counts per sample
        top_clone_marks: `.head` arguments of `repClonoality()`
        top_clone_by: Groupings when visualize top clones
        rare_clone_marks: `.bound` arguments of `repClonoality()`
        rare_clone_by: Groupings when visualize rare clones
        hom_clone_marks: `.clone.types` arguments of `repClonoality()`
        hom_clone_by: Groupings when visualize homeo clones
        overlap_methods: The methods used for `repOverlap()`, each will
            generate a heatmap.
        overlap_redim: Plot the samples with these dimension reduction methods
        gu_by: Groupings to show gene usages
            Multiple groups supported, for example:
            `volume_by = {0: "Status", 1: ["Status", "Sex"]}`
            Or label the groups:
            `volume_by = {"Status": "Status", "Status_Sex": ["Status","Sex"]}`
            If a list or a single variable is given, it will be changed
            into `{"Status": "Status"}`
        gu_top: How many top (ranked by total usage across samples) genes to
            show in the plots
        gua_methods: controls how the data is going to be preprocessed and
            analysed. One of js, cor, cosine, pca, mds, and tsne
        spect: `.quant` and `.col` for `spectratype()` for each sample
        div_methods: Methods to calculate diversities
        div_by: Groupings to show sample diversities
        raref_by: Groupings to show rarefactions
        tracking_target: and
        tracking_samples: The target and samples to track.
            You can do multiple trackings. To do that, you need to specify
            a key for each tracking. It will use the target and samples under
            the same key. If samples from `tracking_samples` cannot be found,
            all samples will be used
            Other than the target supported by immunarch, you can also specify
            top shared clones. For example:
            `tracking_target = { "top_4": {"TOP": 4} }`
        kmers: Arguments for kmer analysis.
            Keys are the K of mers. Values are parameters:
            - `head` specifies # of the most abundant kmers to visualise.
            - `position`: positions of bars: `stack`, `dodge` and `fill`
            - `log`: log-transformation of y-axis
            - `motif`: Method for motif analysis
            There can be multiple `head`s and `motif`s.
            If you do want multiple parameter sets for the same K, You can use
            a float number as the K. For example: `5.1` for K `5`.
    """

    input = "immdata:file"
    output = "outdir:dir:{{in.immdata | stem}}.immunarch"
    lang = config.lang.rscript
    envs = {
        # basic statistics
        "volume_by": {},
        "len_by": {},
        "count_by": {},
        # clonality
        "top_clone_marks": [10, 100, 1000, 3000, 10000],
        "top_clone_by": {},
        "rare_clone_marks": [1, 3, 10, 30, 100],
        "rare_clone_by": {},
        "hom_clone_marks": dict(
            Small=0.0001,
            Medium=0.001,
            Large=0.01,
            Hyperexpanded=1,
        ),
        "hom_clone_by": {},
        # overlapping
        "overlap_methods": ["public"],
        "overlap_redim": ["tsne", "mds"],
        # gene usage
        "gu_by": {},
        "gu_top": 30,
        # gene usage analysis
        "gua_methods": ["js", "cor"],
        # Spectratyping
        "spect": [dict(quant="id", col="nt"), dict(quant="count", col="aa+v")],
        # Diversity
        "div_methods": ["div", "gini.simp"],
        "div_by": {},
        "raref_by": {},
        # Clonotype tracking
        "tracking_target": {},
        "tracking_samples": {},  # can specify order
        # Kmer analysis
        "kmers": {
            5: {"head": 10, "position": "stack", "log": False, "motif": "self"}
        },
    }
    script = "file://../scripts/tcr/Immunarch.R"
    plugin_opts = {
        "report": "file://../reports/tcr/Immunarch.svelte",
        "report_paging": 3,
    }


class CloneResidency(Proc):
    """Identification of clone residency

    Typically, where the clones are located for the sample patient.

    Input:
        immdata: The data loaded by `immunarch::repLoad()`

    Output:
        outdir: The output directory

    Envs:
        subject: The key of subject in metadata. The clone residency will
            be examined for this subject/patient
        group: The key of group in metadata. This usually marks the samples
            that you want to compare. For example, Tumor vs Normal,
            post-treatment vs baseline
            It doesn't have to be 2 groups always. If there are more than 3
            groups, instead of venn diagram, upset plots will be used.
        order: The order of the values in `group`. Early-ordered group will
            be used as x-axis in scatter plots
            If there are more than 2 groups, for example, [A, B, C], the
            scatter plots will be drawn for pairs: B ~ A, C ~ B.
    """

    input = "immdata:file"
    output = "outdir:dir:{{in.immdata | stem}}.cloneov"
    lang = config.lang.rscript
    envs = {
        "subject": [],
        "group": None,
        "order": [],
    }
    script = "file://../scripts/tcr/CloneResidency.R"
    order = 2
    plugin_opts = {"report": "file://../reports/tcr/CloneResidency.svelte"}


class Immunarch2VDJtools(Proc):
    """Convert immuarch format into VDJtools input formats

    Input:
        immdata: The data loaded by `immunarch::repLoad()`

    Output:
        outdir: The output directory containing the vdjtools input for each
            sample
    """
    input = "immdata:file"
    output = "outdir:dir:{{in.immdata | stem}}.vdjtools_input"
    lang = config.lang.rscript
    script = "file://../scripts/tcr/Immunarch2VDJtools.R"


class VJUsage(Proc):
    """Circos-style V-J usage plot displaying the frequency of
    various V-J junctions.

    Input:
        infile: The input file, in vdjtools input format

    Output:
        outfile: The V-J usage plot

    Envs:
        vdjtools: The path to vdjtools
        vdjtools_patch: A patch for vdjtools
    """
    input = "infile:file"
    output = "outfile:file:{{ in.infile | stem | replace: '.vdjtools', '' }}.fancyvj.wt.png"
    lang = config.lang.rscript
    envs = {
        "vdjtools": config.exe.vdjtools,
        "vdjtools_patch": str(SCRIPT_DIR / "tcr" / "vdjtools-patch.sh"),
    }
    order = 3
    script = "file://../scripts/tcr/VJUsage.R"
    plugin_opts = {"report": "file://../reports/tcr/VJUsage.svelte"}


class Attach2Seurat(Proc):
    """Attach the clonal information to a Seurat object as metadata

    Input:
        immfile: The immunarch object in RDS
        sobjfile: The Seurat object file in RDS

    Output:
        outfile: The Seurat object with the clonal information as metadata

    Envs:
        prefix: The prefix to the barcodes. You can use placeholder like
            `{Sample}_` to use the meta data from the immunarch object
        metacols: Which meta columns to attach
    """
    input = "immfile:file, sobjfile:file"
    output = "outfile:file:{{in.sobjfile | basename}}"
    lang = config.lang.rscript
    envs = {
        "prefix": "{Sample}_",
        "metacols": ["Clones", "Proportion", "CDR3.aa"],
    }
    script = "file://../scripts/tcr/Attach2Seurat.R"


class TCRClustering(Proc):
    """Cluster the TCR clones by their CDR3 sequences

    With GIANA

    https://github.com/s175573/GIANA

    > Zhang, Hongyi, Xiaowei Zhan, and Bo Li.
    > "GIANA allows computationally-efficient TCR clustering and multi-disease
    > repertoire classification by isometric transformation."
    > Nature communications 12.1 (2021): 1-11.

    Or ClusTCR

    https://github.com/svalkiers/clusTCR

    > Sebastiaan Valkiers, Max Van Houcke, Kris Laukens, Pieter Meysman,
    > ClusTCR: a Python interface for rapid clustering of large sets of CDR3
    > sequences with unknown antigen specificity,
    > Bioinformatics, 2021.

    Input:
        immfile: The immunarch object in RDS

    Output:
        immfile: The immnuarch object in RDS with TCR cluster information
        clusterfile: The cluster file.
            Columns are CDR3.aa, TCR_Cluster

    Envs:
        tool: The tool used to do the clustering, either GIANA or ClusTCR
            For GIANA, using TRBV mutations is not supported
        on_multi: Whether to run clustering on multi-chain seq or
            the seq read and processed by immunarch
        python: The path of python with `GIANA`'s dependencies installed
            or with `clusTCR` installed. Depending on the `tool` you choose.
        tmpdir: The temporary directory to store the GIANA sources
        giana_repo: The URL prefix for the source code of GIANA
        args: The arguments for the clustering tool
            For GIANA, they will be passed to `python GIAna.py`
            For ClusTCR, they will be passed to `clustcr.Clustering(...)`

    Requires:
        - name: clusTCR
          if: {{ proc.envs.tool == 'ClusTCR' }}
          check: |
            {{ proc.envs.python }} -c "import clustcr"
    """
    input = "immfile:file"
    output = [
        "immfile:file:{{in.immfile | basename}}",
        "clusterfile:file:{{in.immfile | stem}}.clusters.txt",
    ]
    lang = config.lang.rscript
    envs = {
        "tool": "GIANA",  # or ClusTCR
        "on_multi": False,
        "python": config.lang.python,
        "tmpdir": config.path.tmpdir,
        "giana_repo": (
            "https://raw.githubusercontent.com/s175573/GIANA/master/"
        ),
        "args": {},
    }
    script = "file://../scripts/tcr/TCRClustering.R"


class TCRClusteringStats(Proc):
    """Statistics of TCR clusters, generated by biopipen.ns.tcr.TCRClustering

    Input:
        immfile: The immunarch object with TCR clusters attached

    Output:
        outdir: The output directory containing the stats and reports

    Envs:
        shared_clusters: Stats about shared TCR clusters
            numbers_on_heatmap: Whether to show the numbers on the heatmap
            heatmap_meta: The metadata to show on the heatmap
            grouping: The groups to investigate the shared clusters
        sample_diversity: Sample diversity using TCR clusters instead of clones
            keys are the methods and values, currently, `by` to plot
            the diversities by groups

    Requires:
        - name: r-immunarch
          check: |
            {{proc.lang}} -e "library(immunarch)"
    """
    input = "immfile:file"
    output = "outdir:dir:{{in.immfile | stem}}.tcrclusters_stats"
    lang = config.lang.rscript
    envs = {
        "shared_clusters": {
            "numbers_on_heatmap": True,
            "heatmap_meta": [],
            "grouping": None,
        },
        "sample_diversity": {
            "gini": {},  # by = ["Status", "Sex"]
            "gini.simp": {},
        },
    }
    script = "file://../scripts/tcr/TCRClusteringStats.R"
    plugin_opts = {
        "report": "file://../reports/tcr/TCRClusteringStats.svelte",
    }
