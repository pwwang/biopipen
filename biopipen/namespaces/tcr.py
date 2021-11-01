"""Tools to analyze single-cell TCR sequencing data"""
from pipen import Pipen
from pipen.channel import expand_dir

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
        allclones: A text file with all clones across all samples with
            the barcodes for each sample

    Envs:
        tmpdir: The temporary directory to link all data files.
    """
    input = "metafile:file"
    output = [
        "rdsfile:file:{{in.metafile | stem}}.immunarch.RDS",
        "allclones:file:{{in.metafile | stem}}.allclones.txt"
    ]
    lang = config.lang.rscript
    envs = { "tmpdir": config.path.tmpdir }
    script = "file://../scripts/tcr/ImmunarchLoading.R"


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
            Small=.0001,
            Medium=.001,
            Large=.01,
            Hyperexpanded=1,
        ),
        "hom_clone_by": {},
        # overlapping
        "overlap_methods": ["public"],
        "overlap_redim": ["tsne", "mds"],
        # # gene usage
        # "gu_by": [],
        # # gene usage analysis
        # "gua_methods": [],
        # # Spectratyping
        # "spect": [dict(quant="id", col="nt"), dict(quant="count", col="aa+v")],
        # # Diversity
        # "div_methods": [],
        # "div_by": {},
        # "raref_by": [],
        # # Clonotype tracking
        # "tracking_target": [],
        # "tracking_samples": [], # can specify order
        # # Kmer analysis
        # "kmer_len": [5],
        # "kmer_heads": {5: [10, 20]},
        # "kmer_positions": {5: ["dodge", "fill"]},
        # "kmer_log": {}, # False by default
        # # motif analysis
        # "motif_methods": ["self"],
    }
    script = "file://../scripts/tcr/ImmunarchBasic.R"
    plugin_opts = {
        "report": "file://../reports/tcr/ImmunarchBasic.svelte"
    }


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
            `volume_by = {0: "Status", 1: ["Status", "Sex"]}`
            Or label the groups:
            `volume_by = {"Status": "Status", "Status_Sex": ["Status", "Sex"]}`
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
        tracking_samples: The target and samples to tracking
            You can do multiple trackings. To do that, you need to specify
            a key for each tracking. It will use the target and samples under
            the same key. If samples from `tracking_samples` cannot be found,
            all samples will be used
            Other than the target supported by immunarch, you can also specify
            top shared clones. For example:
            `tracking_target = {"top_4": {"TOP": 4}}`
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
        "tracking_samples": {}, # can specify order
        # Kmer analysis
        "kmer_len": [5],
        "kmer_heads": {5: [10, 20]},
        "kmer_positions": {5: ["dodge", "fill"]},
        "kmer_log": {}, # False by default
        # motif analysis
        "motif_methods": ["self"],
    }
    script = "file://../scripts/tcr/ImmunarchAdvanced.R"
    order = 1
    plugin_opts = {
        "report": "file://../reports/tcr/ImmunarchAdvanced.svelte"
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
    plugin_opts = {
        "report": "file://../reports/tcr/CloneResidency.svelte"
    }
