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

    Envs:
        tmpdir: The temporary directory to link all data files.
    """

    input = "metafile:file"
    output = "rdsfile:file:{{in.metafile | stem}}.immunarch.RDS"
    lang = config.lang.rscript
    envs = {"tmpdir": config.path.tmpdir}
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
        groupfile: Also a group file with first column the groups and other
            columns the cell barcodes in the samples

    Envs:
        merge: Merge the cells from the samples, instead of list cells
            for different samples. The cell ids will be preficed with the sample
            name, connected with `_`. The column name will be `ALL` instead.
        clonotype: Use clonotype (CDR3.aa) as the group. The name from
            `envs.filters` will be ignored
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
            For example:
            >>> {{
            >>>   "Top20BM_Post": {{
            >>>     "by.meta": {{"Source": "BM", "Status": "Post"}},
            >>>     "by.count": {{
            >>>         "ORDER": 1, "filter": "TOTAL %in% TOTAL[1:20]"
            >>>     }}
            >>>   }}
            >>> }}
    """
    input = "immdata:file, filterfile:file"
    output = [
        "outfile:file:{{in.immdata | stem}}.RDS",
        "groupfile:file:{{in.immdata | stem}}.groups.txt"
    ]
    envs = {
        "merge": False,
        "clonotype": False,
        "filters": {},
    }
    lang = config.lang.rscript
    script = "file://../scripts/tcr/ImmunarchFilter.R"


class ImmunarchBasic(Proc):
    """Immunarch - Basic statistics and clonality

    See https://immunarch.com/articles/web_only/v3_basic_analysis.html

    Input:
        immdata: The data loaded by `immunarch::repLoad()`

    Output:
        outdir: The output directory

    Envs:
        volume_by: Groupings to show clonotype volume (sizes)
            Multiple groups supported, for example:
            `volume_by = {{0: "Status", 1: ["Status", "Sex"]}}`
            Or label the groups:
            `volume_by = {{"Status": "Status", "Status_Sex": ["Status", "Sex"]}}`
            If a list or a single variable is given, it will be changed
            into `{{"Status": "Status"}}`
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
    """

    input = "immdata:file"
    output = "outdir:dir:{{in.immdata | stem}}.basic"
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
    }
    script = "file://../scripts/tcr/ImmunarchBasic.R"
    plugin_opts = {"report": "file://../reports/tcr/ImmunarchBasic.svelte"}


class ImmunarchAdvanced(Proc):
    """Immunarch - Advanced analysis

    Including gene usage, diversity estimation, tracking clonotype changes and
    kmer/motif analysis

    Input:
        immdata: The data loaded by `immunarch::repLoad()`

    Output:
        outdir: The output directory

    Envs:
        gu_by: Groupings to show gene usages
            Multiple groups supported, for example:
            `volume_by = {{0: "Status", 1: ["Status", "Sex"]}}`
            Or label the groups:
            `volume_by = {{"Status": "Status", "Status_Sex": ["Status","Sex"]}}`
            If a list or a single variable is given, it will be changed
            into `{{"Status": "Status"}}`
        gu_top: How many top (ranked by total usage across samples) genes to
            show in the plots
        gua_methods: controls how the data is going to be preprocessed and
            analysed. One of js, cor, cosine, pca, mds, and tsne
        spect: `.quant` and `.col` for `spectratype()` for each sample
        div_methods: Methods to calculate diversities
        div_by: Groupings to show sample diversities
        raref_by: Groupings to show rarefactions
        tracking_target: and
        tracking_samples: The target and samples to tracking
            You can do multiple trackings. To do that, you need to specify
            a key for each tracking. It will use the target and samples under
            the same key. If samples from `tracking_samples` cannot be found,
            all samples will be used
            Other than the target supported by immunarch, you can also specify
            top shared clones. For example:
            `tracking_target = {{ "top_4": {{"TOP": 4}} }}`
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
    output = "outdir:dir:{{in.immdata | stem}}.advanced"
    lang = config.lang.rscript
    envs = {
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
    script = "file://../scripts/tcr/ImmunarchAdvanced.R"
    order = 1
    plugin_opts = {"report": "file://../reports/tcr/ImmunarchAdvanced.svelte"}


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
