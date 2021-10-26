"""Tools to analyze single-cell TCR sequencing data"""
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
    """
    input = "metafile:file"
    output = "rdsfile:file:{{in.metafile | stem}}.immunarch.RDS"
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
    output = "outdir:dir:{{in.immdata | stem}}.clonality"
    lang = config.lang.rscript
    envs = {
        # basic statistics
        "volume_by": {},
        "len_by": {},
        "count_by": {},
        # clonality
        "top_clone_marks": [10, 100, 1000, 3000, 10000],
        "top_clone_by": {},
        "rare_clone_marks": [1, 3, 10, 30, 100, None],
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
