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
            * `Sample` to specify the sample names.
            * `TCRDir` to assign the path of the data to the samples,
            and this column will be excluded as metadata.
            Immunarch is able to fetch the sample names from the names of
            the target files. However, 10x data yields result like
            `filtered_contig_annotations.csv`, which doesn't have any name
            information.

    Output:
        rdsfile: The RDS file with the data and metadata
        metatxt: The meta data of the cells, used to attach to the Seurat object

    Envs:
        prefix: The prefix to the barcodes. You can use placeholder like
            `{Sample}_` to use the meta data from the immunarch object
        tmpdir (hidden): The temporary directory to link all data files.
        mode (hidden): Either "single" for single chain data or "paired" for
            paired chain data. For `single`, only TRB chain will be kept
            at `immdata$data`, information for other chains will be
            saved at `immdata$tra` and `immdata$multi`.
        metacols (type=list; hidden): The columns to be exported from the
            metatxt.
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

    Analyses include -
    - basic statistics, provided by [`immunarch::repExplore`](https://immunarch.com/reference/repExplore.html)
    such as number of clones or distributions of lengths and counts.
    - the clonality of repertoires, provided by [`immunarch::repClonality`](https://immunarch.com/reference/repClonality.html)
    - the repertoire overlap, provided by [`immunarch::repOverlap`](https://immunarch.com/reference/repOverlap.html)
    - the repertoire overlap, including different clustering procedures and PCA, provided by [`immunarch::repOverlapAnalysis`](https://immunarch.com/reference/repOverlapAnalysis.html)
    - the distributions of V or J genes, provided by [`immunarch::geneUsage`](https://immunarch.com/reference/geneUsage.html)
    - the diversity of repertoires, provided by [`immunarch::repDiversity`](https://immunarch.com/reference/repDiversity.html)
    - the dynamics of repertoires across time points/samples, provided by [`immunarch::trackClonotypes`](https://immunarch.com/reference/trackClonotypes.html)
    - the spectratype of clonotypes, provided by [`immunarch::spectratype`](https://immunarch.com/reference/spectratype.html)
    - the distributions of kmers and sequence profiles, provided by [`immunarch::getKmers`](https://immunarch.com/reference/getKmers.html)

    Input:
        immdata: The data loaded by `immunarch::repLoad()`

    Output:
        outdir: The output directory

    Envs:
        volume_by (ctype=auto): Groupings to show clonotype volume (sizes)
            >>> exp_vol <- repExplore(immdata$data, .method = "volume")
            >>> vis(exp_vol, .by = c("Status"), .meta = immdata$meta)
            >>> vis(exp_vol, .by = c("Status", "Sex"), .meta = immdata$meta)
            This argument will be passed to `.by` of `vis(exp_vol, ...)`.
            It will be finally compiled into a dict with keys the names of
            cases and values the groupings, such as `"Status"` or
            `["Status", "Sex"]`. The value of this argument can be:
            a dict, for example, `{"By_Status": "Status"}`
            a list of strings, for example, `["Status", "Sex"]`, which will
            be compiled to `{"Status_Sex": ["Status", "Sex"]}`
            a string, for example, `"Status"`, which will be compiled to
            `{"Status": "Status"}`
        len_by (ctype=auto): Groupings to show CDR3 length of both aa and nt
            Supported types of values are the same as `volume_by`.
        count_by (ctype=auto): Groupings to show clonotype counts per sample
            Supported types of values are the same as `volume_by`.
        top_clone_marks (ctype=list): A numerical vector with ranges
            of the top clonotypes. Passed to the `.head` argument of
            `repClonoality()`
        top_clone_by (type=auto): Groupings when visualize top clones.
            Supported types of values are the same as `volume_by`.
        rare_clone_marks (ctype=list): A numerical vector with ranges of
            abundance for the rare clonotypes in the dataset. Passed to
            the `.bound` argument of `repClonoality()`
        rare_clone_by (ctype=auto): Groupings when visualize rare clones
            Supported types of values are the same as `volume_by`.
        hom_clone_marks (ctype=ns): A dict with the threshold of the
            half-closed intervals that mark off clonal groups. Passed to the
            `.clone.types` arguments of `repClonoality()`. The keys could be:
            - Rare: the rare clonotypes
            - Small: the small clonotypes
            - Medium: the medium clonotypes
            - Large: the large clonotypes
            - Hyperexpanded: the hyperexpanded clonotypes
        hom_clone_by (ctype=auto): Groupings when visualize homeo clones
            Supported types of values are the same as `volume_by`.
        overlap_methods (ctype=mchoice): The methods used for `repOverlap()`,
            each will generate a heatmap.
            - public: number of public clonotypes between two samples.
            - overlap: a normalised measure of overlap similarity.
                It is defined as the size of the intersection divided by the
                smaller of the size of the two sets.
            - jaccard: conceptually a percentage of how many objects two sets
                have in common out of how many objects they have total.
            - tversky: an asymmetric similarity measure on sets that compares
                a variant to a prototype.
            - cosine: a measure of similarity between two non-zero vectors of
                an inner product space that measures the cosine of the angle
                between them.
            - morisita: how many times it is more likely to randomly select two
                sampled points from the same quadrat (the dataset is covered by
                a regular grid of changing size) then it would be in the case
                of a random distribution generated from a Poisson process.
                Duplicate objects are merged with their counts are summed up.
            - inc+public: incremental overlaps of the N most abundant clonotypes
                with incrementally growing N using the public method.
            - inc+morisita: incremental overlaps of the N most abundant
                clonotypes with incrementally growing N using the morisita
                method.
        overlap_redim (ctype=list): Plot the samples with these dimension
            reduction methods. The methods could be `hclust`, `tsne` or `mds`.
            They could also be combined, for example, `mds+hclust`.
            See https://immunarch.com/reference/repOverlapAnalysis.html
        gu_by (ctype=auto): Groupings to show gene usages
            Supported types of values are the same as `volume_by`.
        gu_top (type=int): How many top (ranked by total usage across samples)
            genes to show in the plots
        gua_methods (type=list): controls how the data is going to be
            preprocessed and analysed. One of js, cor, cosine, pca, mds,
            and tsne. Can also be combined with following methods for the
            actual analysis: hclust, kmeans, dbscan, kruskal. For example:
            `cosine+hclust`.
            See https://immunarch.com/articles/web_only/v5_gene_usage.html
        spect (ctype=json): A list of values for `.quant` and `.col` for
            `spectratype()` for each sample.
        div_methods (mchoice): Methods to calculate diversities
            - chao1: a nonparameteric asymptotic estimator of species richness
                (number of species in a population).
            - hill: Hill numbers are a mathematically unified family of
                diversity indices (differing only by an exponent q).
            - div: true diversity, or the effective number of types, refers to
                the number of equally abundant types needed for the average
                proportional abundance of the types to equal that observed in
                the dataset of interest where all types may not be equally
                abundant.
            - gini.simp: The Gini-Simpson index is the probability of
                interspecific encounter, i.e., probability that two entities
                represent different types.
            - inv.simp: Inverse Simpson index is the effective number of types
                that is obtained when the weighted arithmetic mean is used to
                quantify average proportional abundance of types in the dataset
                of interest.
            - gini: The Gini coefficient measures the inequality among values
                of a frequency distribution (for example levels of income).
                A Gini coefficient of zero expresses perfect equality,
                where all values are the same (for example, where everyone has
                the same income). A Gini coefficient of one (or 100 percents)
                expresses maximal inequality among values (for example where
                only one person has all the income).
            - raref: a technique to assess species richness from the results of
                sampling through extrapolation.
        div_by (ctype=auto): Groupings to show sample diversities
            Supported types of values are the same as `volume_by`.
        raref (ns): Parameters to control the rarefaction analysis
            - by: The variables to group samples
            - separate_by: The variable to separate samples, which will be
                plotted in separate figures. Currently only support one variable
            - align_y (action=store_true): Align max of y-axis if there are
                multiple figures
            - align_x (action=store_true): Align max of x-axis if there are
                multiple figures
            - log (action=store_true): Also plot log-transformed x-axis using
                `vis(.log = TRUE)`.
            - <other>: Other arguments for `repDiversity(.method="raref", ...)`
                i.e. ".step", ".norm"
        tracking_target (ctype=json): Either a list of AA seq of clonotypes
            to track, or a dict of those lists. The keys will be used as
            the names of the tracks. If you want to track the top N clonotypes,
            you can use `{"TOP": N}`.
        tracking_samples (ctype=json): The samples to track. If not specified,
            all samples will be used. Make sure the keys in `tracking_target`
            and `tracking_samples` are the same, if you want to track multiple
            cases at the same time.
        kmers (ctype=json): Arguments for kmer analysis.
            There can be multiple `head`s and `motif`s.
            If you do want multiple parameter sets for the same K, You can use
            a float number as the K. For example: `5.1` for K `5`.
            Keys are the K of mers. Values are parameters
            - head: specifies # of the most abundant kmers to visualise.
            - position: positions of bars: `stack`, `dodge` and `fill`
            - log: log-transformation of y-axis
            - motif: Method for motif analysis
    """  # noqa: E501
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
        "raref": {
            "by": None,
            "separate_by": None,
            "align_x": False,
            "align_y": False,
            "log": False,
            # ... other arguments for repDiversity()
            # i.e. ".step", ".norm"
        },
        # Clonotype tracking
        "tracking_target": {},
        "tracking_samples": {},  # can specify order
        # Kmer analysis
        "kmers": {
            "5": {"head": 10, "position": "stack", "log": False, "motif": "self"}
        },
    }
    script = "file://../scripts/tcr/Immunarch.R"
    plugin_opts = {
        "report": "file://../reports/tcr/Immunarch.svelte",
        "report_paging": 3,
    }


class SampleDiversity(Proc):
    """Sample diversity and rarefaction analysis

    This is part of Immunarch, in case we have multiple dataset to compare.

    Input:
        immdata: The data loaded by `immunarch::repLoad()`

    Output:
        outdir: The output directory

    Envs:
        div_methods: Methods to calculate diversities
            It is a dict, keys are the method names, values are the groupings.
            Each one is a case, multiple columns for a case are separated by `,`
            For example: `{"div": ["Status", "Sex", "Status,Sex"]}` will run
            true diversity for samples grouped by `Status`, `Sex`, and both.
            The diversity for each sample without grouping will also be added
            anyway.
            Supported methods: `chao1`, `hill`, `div`, `gini.simp`, `inv.simp`,
            `gini`, and `raref`. See also
            https://immunarch.com/articles/web_only/v6_diversity.html
        devpars: The parameters for the plotting device
            It is a dict, and keys are the methods and values are dicts with
            width, height and res that will be passed to `png()`
            If not provided, 1000, 1000 and 100 will be used.
    """
    input = "immdata:file"
    output = "outdir:dir:{{in.immdata | stem}}.diversity"
    lang = config.lang.rscript
    envs = {
        "div_methods": {
            "chao1": [],
            "hill": [],
            "div": [],
            "gini.simp": [],
            "inv.simp": [],
            "gini": [],
            "raref": [],
        },
        "devpars": {},
    }
    script = "file://../scripts/tcr/SampleDiversity.R"
    plugin_opts = {
        "report": "file://../reports/tcr/SampleDiversity.svelte",
    }


class CloneResidency(Proc):
    """Identification of clone residency

    Typically, where the clones are located for the sample patient.

    Input:
        immdata: The data loaded by `immunarch::repLoad()`

    Output:
        outdir: The output directory

    Envs:
        subject (ctyle=list): The key of subject in metadata. The clone
            residency will be examined for this subject/patient
        group: The key of group in metadata. This usually marks the samples
            that you want to compare. For example, Tumor vs Normal,
            post-treatment vs baseline
            It doesn't have to be 2 groups always. If there are more than 3
            groups, instead of venn diagram, upset plots will be used.
        order (ctyle=list): The order of the values in `group`. Early-ordered
            group will be used as x-axis in scatter plots
            If there are more than 2 groups, for example, [A, B, C], the
            scatter plots will be drawn for pairs: B ~ A, C ~ B and C ~ A.
        sample_groups: How the samples aligned in the report.
            Useful for cohort with large number of samples.
    """
    input = "immdata:file"
    output = "outdir:dir:{{in.immdata | stem}}.cloneov"
    lang = config.lang.rscript
    envs = {
        "subject": [],
        "group": None,
        "order": [],
        "sample_groups": None,
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


class ImmunarchSplitIdents(Proc):
    """Split the data into multiple immunarch datasets by Idents from Seurat

    Note that only the cells in both the `immdata` and `sobjfile` will be
    kept.

    Requires `immunarch >= 0.9.0` to use `select_clusters()`

    Input:
        immdata: The data loaded by `immunarch::repLoad()`
        sobjfile: The Seurat object file.
            You can set a different ident by `Idents(sobj) <- "new_ident"` to
            split the data by the new ident, where `"new_ident"` is the an
            existing column in meta data

    Output:
        outdir: The output directory containing the RDS files of the splitted
            immunarch datasets

    Envs:
        prefix: The prefix of the cell barcodes in the `Seurat` object.
            Once could use a fixed prefix, or a placeholder with the column
            name in meta data. For example, `"{Sample}_"` will replace the
            placeholder with the value of the column `Sample` in meta data.
        sample_col: The column name in meta data that contains the sample name
    """
    input = "immdata:file, sobjfile:file"
    output = "outdir:dir:{{in.immdata | stem}}.splitidents"
    lang = config.lang.rscript
    envs = {"prefix": "{Sample}_", "sample_col": "Sample"}
    script = "file://../scripts/tcr/ImmunarchSplitIdents.R"


class VJUsage(Proc):
    """Circos-style V-J usage plot displaying the frequency of
    various V-J junctions using vdjtools

    Input:
        infile: The input file, in vdjtools input format

    Output:
        outfile: The V-J usage plot

    Envs:
        vdjtools: The path to vdjtools
        vdjtools_patch (hidden): A patch for vdjtools
    """

    input = "infile:file"
    output = (
        "outfile:file:{{ in.infile | stem | replace: '.vdjtools', '' }}"
        ".fancyvj.wt.png"
    )
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
        tool (choice): The tool used to do the clustering, either
            [GIANA](https://github.com/s175573/GIANA) or
            [ClusTCR](https://github.com/svalkiers/clusTCR).
            For GIANA, using TRBV mutations is not supported
            - GIANA: by Li lab at UT Southwestern Medical Center
            - ClusTCR: by Sebastiaan Valkiers, etc
        on_multi (action=store_true;hidden): Whether to run clustering on
            multi-chain seq or the seq read and processed by immunarch
        python: The path of python with `GIANA`'s dependencies installed
            or with `clusTCR` installed. Depending on the `tool` you choose.
        tmpdir: The temporary directory to store the GIANA sources
        giana_repo: The URL prefix for the source code of GIANA
        args (type=json): The arguments for the clustering tool
            For GIANA, they will be passed to `python GIAna.py`
            See https://github.com/s175573/GIANA#usage
            For ClusTCR, they will be passed to `clustcr.Clustering(...)`
            See https://svalkiers.github.io/clusTCR/docs/clustering/how-to-use.html#clustering

    Requires:
        clusTCR:
            - if: {{ proc.envs.tool == 'ClusTCR' }}
            - check: {{ proc.envs.python }} -c "import clustcr"
    """  # noqa: E501
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
        shared_clusters (ns): Stats about shared TCR clusters
            - numbers_on_heatmap (action=store_true): Whether to show the
                numbers on the heatmap
            - heatmap_meta (list): The metadata to show on the heatmap
            - grouping: The groups to investigate the shared clusters
        sample_diversity (ctype=json): Sample diversity using TCR clusters
            instead of clones keys are the methods and values, currently, `by`
            to plot the diversities by groups

    Requires:
        r-immunarch:
            - check: {{proc.lang}} -e "library(immunarch)"
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


class CloneSizeQQPlot(Proc):
    """QQ plot of the clone sizes

    QQ plots for clones sizes of pairs of samples

    Input:
        immdata: The data loaded by `immunarch::repLoad()`

    Output:
        outdir: The output directory

    Envs:
        subject: The key of subject in metadata, defining the pairs.
            The clone residency will be examined for this subject/patient
        group: The key of group in metadata. This usually marks the samples
            that you want to compare. For example, Tumor vs Normal,
            post-treatment vs baseline
            It doesn't have to be 2 groups always. If there are more than 3
            groups, for example, [A, B, C], the QQ plots will be generated
            for all the combinations of 2 groups, i.e., [A, B], [A, C], [B, C]
        order: The order of the values in `group`. Early-ordered group will
            be used as x-axis in scatter plots
            If there are more than 2 groups, for example, [A, B, C], the
            QQ plots will be drawn for pairs: B ~ A, C ~ B.
        diag: Whether to draw the diagonal line in the QQ plot
        on: The key of the metadata to use for the QQ plot. One/Both of
            `["Clones", "Proportion"]`
    """
    input = "immdata:file"
    output = "outdir:dir:{{in.immdata | stem}}.qqplots"
    lang = config.lang.rscript
    envs = {
        "subject": [],
        "group": None,
        "order": [],
        "diag": True,
        "on": ["Clones", "Proportion"],
    }
    script = "file://../scripts/tcr/CloneSizeQQPlot.R"
    order = 3
    plugin_opts = {"report": "file://../reports/tcr/CloneSizeQQPlot.svelte"}
