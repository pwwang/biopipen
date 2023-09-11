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
            * `TCRData` to assign the path of the data to the samples,
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
        metacols (list; hidden): The columns to be exported to the metatxt.
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
        mutaters (type=json;order=-9): The mutaters passed to `dplyr::mutate()` on `immdata$meta` to add new columns.
            The keys will be the names of the columns, and the values will be the expressions.
            The new names can be used in `volumes`, `lens`, `counts`, `top_clones`, `rare_clones`, `hom_clones`, `gene_usages`, `divs`, etc.
        volumes (ns): Explore clonotype volume (sizes).
            - by: Groupings when visualize clonotype volumes, passed to the `.by` argument of `vis(imm_vol, .by = <values>)`.
                Multiple columns should be separated by `,`.
            - devpars (ns): The parameters for the plotting device.
                - width (type=int): The width of the plot.
                - height (type=int): The height of the plot.
                - res (type=int): The resolution of the plot.
            - cases (type=json;order=9): If you have multiple cases, you can use this argument to specify them.
                The keys will be the names of the cases.
                The values will be passed to the corresponding arguments above.
                If any of these arguments are not specified, the values in `envs.volumes` will be used.
                If NO cases are specified, the default case will be added, with the name `DEFAULT` and the
                values of `envs.volume.by`, `envs.volume.devpars`.
        lens (ns): Explore clonotype CDR3 lengths.
            - by: Groupings when visualize clonotype lengths, passed to the `.by` argument of `vis(imm_len, .by = <values>)`.
                Multiple columns should be separated by `,`.
            - devpars (ns): The parameters for the plotting device.
                - width (type=int): The width of the plot.
                - height (type=int): The height of the plot.
                - res (type=int): The resolution of the plot.
            - cases (type=json;order=9): If you have multiple cases, you can use this argument to specify them.
                The keys will be the names of the cases.
                The values will be passed to the corresponding arguments above.
                If any of these arguments are not specified, the values in `envs.lens` will be used.
                If NO cases are specified, the default case will be added, with the name `DEFAULT` and the
                values of `envs.lens.by`, `envs.lens.devpars`.
        counts (ns): Explore clonotype counts.
            - by: Groupings when visualize clonotype counts, passed to the `.by` argument of `vis(imm_count, .by = <values>)`.
                Multiple columns should be separated by `,`.
            - devpars (ns): The parameters for the plotting device.
                - width (type=int): The width of the plot.
                - height (type=int): The height of the plot.
                - res (type=int): The resolution of the plot.
            - cases (type=json;order=9): If you have multiple cases, you can use this argument to specify them.
                The keys will be the names of the cases.
                The values will be passed to the corresponding arguments above.
                If any of these arguments are not specified, the values in `envs.counts` will be used.
                If NO cases are specified, the default case will be added, with the name `DEFAULT` and the
                values of `envs.counts.by`, `envs.counts.devpars`.
        top_clones (ns): Explore top clonotypes.
            - by: Groupings when visualize top clones, passed to the `.by` argument of `vis(imm_top, .by = <values>)`.
                Multiple columns should be separated by `,`.
            - marks (list;itype=int): A numerical vector with ranges of the top clonotypes. Passed to the `.head` argument of `repClonoality()`.
            - devpars (ns): The parameters for the plotting device.
                - width (type=int): The width of the plot.
                - height (type=int): The height of the plot.
                - res (type=int): The resolution of the plot.
            - cases (type=json;order=9): If you have multiple cases, you can use this argument to specify them.
                The keys will be the names of the cases.
                The values will be passed to the corresponding arguments above.
                If any of these arguments are not specified, the values in `envs.top_clones` will be used.
                If NO cases are specified, the default case will be added, with the name `DEFAULT` and the
                values of `envs.top_clones.by`, `envs.top_clones.marks` and `envs.top_clones.devpars`.
        rare_clones (ns): Explore rare clonotypes.
            - by: Groupings when visualize rare clones, passed to the `.by` argument of `vis(imm_rare, .by = <values>)`.
                Multiple columns should be separated by `,`.
            - marks (list;itype=int): A numerical vector with ranges of abundance for the rare clonotypes in the dataset.
                Passed to the `.bound` argument of `repClonoality()`.
            - devpars (ns): The parameters for the plotting device.
                - width (type=int): The width of the plot.
                - height (type=int): The height of the plot.
                - res (type=int): The resolution of the plot.
            - cases (type=json;order=9): If you have multiple cases, you can use this argument to specify them.
                The keys will be the names of the cases.
                The values will be passed to the corresponding arguments above.
                If any of these arguments are not specified, the values in `envs.rare_clones` will be used.
                If NO cases are specified, the default case will be added, with the name `DEFAULT` and the
                values of `envs.rare_clones.by`, `envs.rare_clones.marks` and `envs.rare_clones.devpars`.
        hom_clones (ns): Explore homeo clonotypes.
            - by: Groupings when visualize homeo clones, passed to the `.by` argument of `vis(imm_hom, .by = <values>)`.
                Multiple columns should be separated by `,`.
            - marks (ns): A dict with the threshold of the half-closed intervals that mark off clonal groups.
                Passed to the `.clone.types` arguments of `repClonoality()`.
                The keys could be:
                - Rare (type=float): the rare clonotypes
                - Small (type=float): the small clonotypes
                - Medium (type=float): the medium clonotypes
                - Large (type=float): the large clonotypes
                - Hyperexpanded (type=float): the hyperexpanded clonotypes
            - devpars (ns): The parameters for the plotting device.
                - width (type=int): The width of the plot.
                - height (type=int): The height of the plot.
                - res (type=int): The resolution of the plot.
            - cases (type=json;order=9): If you have multiple cases, you can use this argument to specify them.
                The keys will be the names of the cases.
                The values will be passed to the corresponding arguments above.
                If any of these arguments are not specified, the values in `envs.hom_clones` will be used.
                If NO cases are specified, the default case will be added, with the name `DEFAULT` and the
                values of `envs.hom_clones.by`, `envs.hom_clones.marks` and `envs.hom_clones.devpars`.
        overlaps (ns): Explore clonotype overlaps.
            - method (choice): The method to calculate overlaps.
                - public: number of public clonotypes between two samples.
                - overlap: a normalised measure of overlap similarity.
                    It is defined as the size of the intersection divided by the smaller of the size of the two sets.
                - jaccard: conceptually a percentage of how many objects two sets have in common out of how many objects they have total.
                - tversky: an asymmetric similarity measure on sets that compares a variant to a prototype.
                - cosine: a measure of similarity between two non-zero vectors of an inner product space that measures the cosine of the angle between them.
                - morisita: how many times it is more likely to randomly select two sampled points from the same quadrat (the dataset is
                    covered by a regular grid of changing size) then it would be in the case of a random distribution generated from
                    a Poisson process. Duplicate objects are merged with their counts are summed up.
                - inc+public: incremental overlaps of the N most abundant clonotypes with incrementally growing N using the public method.
                - inc+morisita: incremental overlaps of the N most abundant clonotypes with incrementally growing N using the morisita method.
            - vis_args (type=json): Other arguments for the plotting functions `vis(imm_ov, ...)`.
            - devpars (ns): The parameters for the plotting device.
                - width (type=int): The width of the plot.
                - height (type=int): The height of the plot.
                - res (type=int): The resolution of the plot.
            - analyses (ns;order=8): Perform overlap analyses.
                - method: Plot the samples with these dimension reduction methods.
                    The methods could be `hclust`, `tsne` or `mds`.
                    They could also be combined, for example, `mds+hclust`.
                    See https://immunarch.com/reference/repOverlapAnalysis.html
                - vis_args (type=json): Other arguments for the plotting functions.
                - devpars (ns): The parameters for the plotting device.
                    - width (type=int): The width of the plot.
                    - height (type=int): The height of the plot.
                    - res (type=int): The resolution of the plot.
                - cases (type=json): If you have multiple cases, you can use this argument to specify them.
                    The keys will be the names of the cases.
                    The values will be passed to the corresponding arguments above.
                    If any of these arguments are not specified, the values in `envs.overlaps.analyses` will be used.
                    If NO cases are specified, the default case will be added, with the name `DEFAULT` and the
                    values of `envs.overlaps.analyses.method`, `envs.overlaps.analyses.vis_args` and `envs.overlaps.analyses.devpars`.
            - cases (type=json;order=9): If you have multiple cases, you can use this argument to specify them.
                The keys will be the names of the cases.
                The values will be passed to the corresponding arguments above.
                If any of these arguments are not specified, the values in `envs.overlaps` will be used.
                If NO cases are specified, the default case will be added, with the key the default method and the
                values of `envs.overlaps.method`, `envs.overlaps.vis_args`, `envs.overlaps.devpars` and `envs.overlaps.analyses`.
        gene_usages (ns): Explore gene usages.
            - top (type=int): How many top (ranked by total usage across samples) genes to show in the plots.
                Use `0` to use all genes.
            - norm (flag): If True then use proportions of genes, else use counts of genes.
            - by: Groupings to show gene usages, passed to the `.by` argument of `vis(imm_gu_top, .by = <values>)`.
                Multiple columns should be separated by `,`.
            - vis_args (type=json): Other arguments for the plotting functions.
            - devpars (ns): The parameters for the plotting device.
                - width (type=int): The width of the plot.
                - height (type=int): The height of the plot.
                - res (type=int): The resolution of the plot.
            - analyses (ns;order=8): Perform gene usage analyses.
                - method: The method to control how the data is going to be preprocessed and analysed.
                    One of `js`, `cor`, `cosine`, `pca`, `mds` and `tsne`. Can also be combined with following methods
                    for the actual analyses: `hclust`, `kmeans`, `dbscan`, and `kruskal`. For example: `cosine+hclust`.
                    See https://immunarch.com/articles/web_only/v5_gene_usage.html.
                - vis_args (type=json): Other arguments for the plotting functions.
                - devpars (ns): The parameters for the plotting device.
                    - width (type=int): The width of the plot.
                    - height (type=int): The height of the plot.
                    - res (type=int): The resolution of the plot.
                - cases (type=json): If you have multiple cases, you can use this argument to specify them.
                    The keys will be the names of the cases.
                    The values will be passed to the corresponding arguments above.
                    If any of these arguments are not specified, the values in `envs.gene_usages.analyses` will be used.
                    If NO cases are specified, the default case will be added, with the name `DEFAULT` and the
                    values of `envs.gene_usages.analyses.method`, `envs.gene_usages.analyses.vis_args` and `envs.gene_usages.analyses.devpars`.
            - cases (type=json;order=9): If you have multiple cases, you can use this argument to specify them.
                The keys will be used as the names of the cases.
                The values will be passed to the corresponding arguments above.
                If any of these arguments are not specified, the values in `envs.gene_usages` will be used.
                If NO cases are specified, the default case will be added, with the name `DEFAULT` and the
                values of `envs.gene_usages.top`, `envs.gene_usages.norm`, `envs.gene_usages.by`, `envs.gene_usages.vis_args`, `envs.gene_usages.devpars` and `envs.gene_usages.analyses`.
        spects (ns): Spectratyping analysis.
            - quant: Select the column with clonal counts to evaluate.
                Set to `id` to count every clonotype once.
                Set to `count` to take into the account number of clones per clonotype.
                Multiple columns should be separated by `,`.
            - col: A string that specifies the column(s) to be processed.
                The output is one of the following strings, separated by the plus sign: "nt" for nucleotide sequences,
                "aa" for amino acid sequences, "v" for V gene segments, "j" for J gene segments.
                E.g., pass "aa+v" for spectratyping on CDR3 amino acid sequences paired with V gene segments,
                i.e., in this case a unique clonotype is a pair of CDR3 amino acid and V gene segment.
                Clonal counts of equal clonotypes will be summed up.
            - devpars (ns): The parameters for the plotting device.
                - width (type=int): The width of the plot.
                - height (type=int): The height of the plot.
                - res (type=int): The resolution of the plot.
            - cases (type=json;order=9): If you have multiple cases, you can use this argument to specify them.
                The keys will be the names of the cases.
                The values will be passed to the corresponding arguments above.
                If any of these arguments are not specified, the values in `envs.spects` will be used.
                By default, a `By_Clonotype` case will be added, with the values of `quant = "id"` and `col = "nt"`, and
                a `By_Num_Clones` case will be added, with the values of `quant = "count"` and `col = "aa+v"`.
        divs (ns): Parameters to control the diversity analysis.
            - filter: The filter passed to `dplyr::filter()` to filter the data for each sample before calculating diversity.
                For example, `Clones > 1` to filter out singletons.
                To check which columns are available, use `immdata$data[[1]] |> colnames()` in R.
                You may also check quickly here:
                https://immunarch.com/articles/v2_data.html#basic-data-manipulations-with-dplyr-and-immunarch
                To use the top 10 clones, you can try `rank(desc(Clones)) <= 10`
            - method (choice): The method to calculate diversity.
                - chao1: a nonparameteric asymptotic estimator of species richness.
                    (number of species in a population).
                - hill: Hill numbers are a mathematically unified family of diversity indices.
                    (differing only by an exponent q).
                - div: true diversity, or the effective number of types.
                    It refers to the number of equally abundant types needed for the average proportional abundance of the types to equal
                    that observed in the dataset of interest where all types may not be equally abundant.
                - gini.simp: The Gini-Simpson index.
                    It is the probability of interspecific encounter, i.e., probability that two entities represent different types.
                - inv.simp: Inverse Simpson index.
                    It is the effective number of types that is obtained when the weighted arithmetic mean is used to quantify
                    average proportional abundance of types in the dataset of interest.
                - gini: The Gini coefficient.
                    It measures the inequality among values of a frequency distribution (for example levels of income).
                    A Gini coefficient of zero expresses perfect equality, where all values are the same (for example, where everyone has the same income).
                    A Gini coefficient of one (or 100 percents) expresses maximal inequality among values (for example where only one person has all the income).
                - d50: The D50 index.
                    It is the number of types that are needed to cover 50%% of the total abundance.
                - dxx: The Dxx index.
                    It is the number of types that are needed to cover xx%% of the total abundance.
                    The percentage should be specified in the `args` argument using `perc` key.
                - raref: Species richness from the results of sampling through extrapolation.
            - by: The variables (column names) to group samples.
                Multiple columns should be separated by `,`.
            - args (type=json): Other arguments for `repDiversity()`.
                Do not include the preceding `.` and use `-` instead of `.` in the argument names.
                For example, `-do-norm` will be compiled to `.do.norm`.
                See all arguments at
                https://immunarch.com/reference/repDiversity.html
            - order (list): The order of the values in `by` on the x-axis of the plots.
                If not specified, the values will be used as-is.
            - test (ns): Perform statistical tests between each pair of groups.
                Does NOT work for `raref`.
                - method (choice): The method to perform the test
                    - none: No test
                    - t.test: Welch's t-test
                    - wilcox.test: Wilcoxon rank sum test
                - padjust (choice): The method to adjust p-values.
                    Defaults to `none`.
                    - bonferroni: one-step correction
                    - holm: step-down method using Bonferroni adjustments
                    - hochberg: step-up method (independent)
                    - hommel: closed method based on Simes tests (non-negative)
                    - BH: Benjamini & Hochberg (non-negative)
                    - BY: Benjamini & Yekutieli (negative)
                    - fdr: Benjamini & Hochberg (non-negative)
                    - none: no correction.
            - separate_by: A column name used to separate the samples into different plots. Only works for `raref`.
            - align_x (flag): Align the x-axis of multiple plots. Only works for `raref`.
            - align_y (flag): Align the y-axis of multiple plots. Only works for `raref`.
            - log (flag): Indicate whether we should plot with log-transformed x-axis using `vis(.log = TRUE)`. Only works for `raref`.
            - devpars (ns): The parameters for the plotting device.
                - width (type=int): The width of the device
                - height (type=int): The height of the device
                - res (type=int): The resolution of the device
            - cases (type=json;order=9): If you have multiple cases, you can use this argument to specify them.
                The keys will be used as the names of the cases.
                The values will be passed to the corresponding arguments above.
                If NO cases are specified, the default case will be added, with the name of `envs.div.method`.
                The values specified in `envs.div` will be used as the defaults for the cases here.
        trackings (ns): Parameters to control the clonotype tracking analysis.
            - targets: Either a set of CDR3AA seq of clonotypes to track (separated by `,`), or simply an integer to track the top N clonotypes.
            - subject_col: The column name in meta data that contains the subjects/samples on the x-axis of the alluvial plot.
                If the values in this column are not unique, the values will be merged with the values in `subject_col` to form the x-axis.
                This defaults to `Sample`.
            - subjects (list): A list of values from `subject_col` to show in the alluvial plot on the x-axis.
                If not specified, all values in `subject_col` will be used.
                This also specifies the order of the x-axis.
            - cases (type=json;order=9): If you have multiple cases, you can use this argument to specify them.
                The keys will be used as the names of the cases.
                The values will be passed to the corresponding arguments (`target`, `subject_col`, and `subjects`).
                If any of these arguments are not specified, the values in `envs.trackings` will be used.
                If NO cases are specified, the default case will be added, with the name `DEFAULT` and the
                values of `envs.trackings.target`, `envs.trackings.subject_col`, and `envs.trackings.subjects`.
        kmers (ns): Arguments for kmer analysis.
            - k (type=int): The length of kmer.
            - head (type=int): The number of top kmers to show.
            - vis_args (type=json): Other arguments for the plotting functions.
            - devpars (ns): The parameters for the plotting device.
                - width (type=int): The width of the plot.
                - height (type=int): The height of the plot.
                - res (type=int): The resolution of the plot.
            - profiles (ns;order=8): Arguments for sequence profilings.
                - method (choice): The method for the position matrix.
                    For more information see https://en.wikipedia.org/wiki/Position_weight_matrix.
                    - freq: position frequency matrix (PFM) - a matrix with occurences of each amino acid in each position.
                    - prob: position probability matrix (PPM) - a matrix with probabilities of each amino acid in each position.
                    - wei: position weight matrix (PWM) - a matrix with log likelihoods of PPM elements.
                    - self: self-information matrix (SIM) - a matrix with self-information of elements in PWM.
                - vis_args (type=json): Other arguments for the plotting functions.
                - devpars (ns): The parameters for the plotting device.
                    - width (type=int): The width of the plot.
                    - height (type=int): The height of the plot.
                    - res (type=int): The resolution of the plot.
                - cases (type=json): If you have multiple cases, you can use this argument to specify them.
                    The keys will be the names of the cases.
                    The values will be passed to the corresponding arguments above.
                    If any of these arguments are not specified, the values in `envs.kmers.profiles` will be used.
                    If NO cases are specified, the default case will be added, with the name `DEFAULT` and the
                    values of `envs.kmers.profiles.method`, `envs.kmers.profiles.vis_args` and `envs.kmers.profiles.devpars`.
            - cases (type=json;order=9): If you have multiple cases, you can use this argument to specify them.
                The keys will be used as the names of the cases.
                The values will be passed to the corresponding arguments above.
                If any of these arguments are not specified, the default case will be added, with the name `DEFAULT` and the
                values of `envs.kmers.k`, `envs.kmers.head`, `envs.kmers.vis_args` and `envs.kmers.devpars`.
    """  # noqa: E501
    input = "immdata:file"
    output = "outdir:dir:{{in.immdata | stem}}.immunarch"
    lang = config.lang.rscript
    envs = {
        "mutaters": {},
        # basic statistics
        "volumes": {
            "by": None,
            "devpars": {"width": 1000, "height": 1000, "res": 100},
            "cases": {},
        },
        "lens": {
            "by": None,
            "devpars": {"width": 1000, "height": 1000, "res": 100},
            "cases": {},
        },
        "counts": {
            "by": None,
            "devpars": {"width": 1000, "height": 1000, "res": 100},
            "cases": {},
        },
        # clonality
        "top_clones": {
            "by": None,
            "marks": [10, 100, 1000, 3000, 10000, 30000, 1e5],
            "devpars": {"width": 1000, "height": 1000, "res": 100},
            "cases": {},
        },
        "rare_clones": {
            "by": None,
            "marks": [1, 3, 10, 30, 100],
            "devpars": {"width": 1000, "height": 1000, "res": 100},
            "cases": {},
        },
        "hom_clones": {
            "by": None,
            "marks": dict(
                Rare=1e-5,
                Small=1e-4,
                Medium=1e-3,
                Large=0.01,
                Hyperexpanded=1.0,
            ),
            "devpars": {"width": 1000, "height": 1000, "res": 100},
            "cases": {},
        },
        # overlapping
        "overlaps": {
            "method": "public",
            "vis_args": {},
            "devpars": {"width": 1000, "height": 1000, "res": 100},
            "analyses": {
                "method": "tsne",
                "vis_args": {},
                "devpars": {"width": 1000, "height": 1000, "res": 100},
                "cases": {},
            },
            "cases": {},
        },
        # gene usage
        "gene_usages": {
            "top": 30,
            "norm": False,
            "by": None,
            "vis_args": {},
            "devpars": {"width": 1000, "height": 1000, "res": 100},
            "analyses": {
                "method": "tsne",
                "vis_args": {},
                "devpars": {"width": 1000, "height": 1000, "res": 100},
                "cases": {},
            },
            "cases": {},
        },
        # Spectratyping
        "spects": {
            "quant": None,
            "col": None,
            "devpars": {"width": 1000, "height": 1000, "res": 100},
            "cases": {
                "By_Clonotype": dict(quant="id", col="nt"),
                "By_Num_Clones": dict(quant="count", col="aa+v"),
            },
        },
        # Diversity
        "divs": {
            "filter": None,
            "method": "gini",
            "by": None,
            "args": {},
            "order": [],
            "test": {
                "method": "none",
                "padjust": "none",
            },
            "separate_by": None,
            "align_x": False,
            "align_y": False,
            "log": False,
            "devpars": {
                "width": 1000,
                "height": 1000,
                "res": 100,
            },
            "cases": {},
        },
        # Clonotype tracking
        "trackings": {
            "targets": None,  # Do not do trackings by default
            "subject_col": "Sample",
            "subjects": [],
            "cases": {},
        },
        # Kmer analysis
        "kmers": {
            "k": 5,
            "head": 10,
            "vis_args": {},
            "devpars": {"width": 1000, "height": 1000, "res": 100},
            "profiles": {
                "method": "self",
                "vis_args": {},
                "devpars": {"width": 1000, "height": 1000, "res": 100},
                "cases": {},
            },
            "cases": {},
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
        subject (list): The key of subject in metadata. The clone
            residency will be examined for this subject/patient
        group: The key of group in metadata. This usually marks the samples
            that you want to compare. For example, Tumor vs Normal,
            post-treatment vs baseline
            It doesn't have to be 2 groups always. If there are more than 3
            groups, instead of venn diagram, upset plots will be used.
        order (list): The order of the values in `group`. Early-ordered
            group will be used as x-axis in scatter plots
            If there are more than 2 groups, for example, [A, B, C], the
            scatter plots will be drawn for pairs: B ~ A, C ~ B and C ~ A.
        sample_groups: How the samples aligned in the report.
            Useful for cohort with large number of samples.
        mutaters (type=json): The mutaters passed to `dplyr::mutate()` on
            `immdata$meta` to add new columns. The keys will be the names of
            the columns, and the values will be the expressions. The new names
            can be used in `subject`, `group`, `order` and `sample_groups`.
        cases (type=json): If you have multiple cases, you can use this argument
            to specify them. The keys will be used as the names of the cases.
            The values will be passed to the corresponding arguments.
            If no cases are specified, the default case will be added, with
            the name `DEFAULT` and the values of `envs.subject`, `envs.group`,
            `envs.order` and `envs.sample_groups`. These values are also the
            defaults for the other cases.
    """
    input = "immdata:file"
    output = "outdir:dir:{{in.immdata | stem}}.cloneov"
    lang = config.lang.rscript
    envs = {
        "subject": [],
        "group": None,
        "order": [],
        "sample_groups": None,
        "mutaters": {},
        "cases": {},
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
        python: The path of python with `GIANA`'s dependencies installed
            or with `clusTCR` installed. Depending on the `tool` you choose.
        args (type=json): The arguments for the clustering tool
            For GIANA, they will be passed to `python GIAna.py`
            See https://github.com/s175573/GIANA#usage
            For ClusTCR, they will be passed to `clustcr.Clustering(...)`
            See https://svalkiers.github.io/clusTCR/docs/clustering/how-to-use.html#clustering
        on_multi (flag;hidden): Whether to run clustering on
            multi-chain seq or the seq read and processed by immunarch

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
        cluster_size (ns): The distribution of size of each cluster.
            - by: The variables (column names) used to fill the histogram.
                Only a single column is supported.
            - devpars (ns): The parameters for the plotting device.
                - width (type=int): The width of the device
                - height (type=int): The height of the device
                - res (type=int): The resolution of the device
            - cases (type=json): If you have multiple cases, you can use this
                argument to specify them. The keys will be the names of the
                cases. The values will be passed to the corresponding arguments
                above. If any of these arguments are not specified, the values
                in `envs.cluster_size` will be used. If NO cases are
                specified, the default case will be added, with the name
                `DEFAULT`.
        shared_clusters (ns): Stats about shared TCR clusters
            - numbers_on_heatmap (flag): Whether to show the
                numbers on the heatmap.
            - heatmap_meta (list): The columns of metadata to show on the
                heatmap.
            - grouping: The groups to investigate the shared clusters.
                If specified, venn diagrams will be drawn instead of heatmaps.
                In such case, `numbers_on_heatmap` and `heatmap_meta` will be
                ignored.
            - devpars (ns): The parameters for the plotting device.
                - width (type=int): The width of the device
                - height (type=int): The height of the device
                - res (type=int): The resolution of the device
            - cases (type=json): If you have multiple cases, you can use this
                argument to specify them. The keys will be the names of the
                cases. The values will be passed to the corresponding arguments
                above. If any of these arguments are not specified, the values
                in `envs.shared_clusters` will be used. If NO cases are
                specified, the default case will be added, with the name
                `DEFAULT`.
        sample_diversity (ns): Sample diversity using TCR clusters instead of
            clones.
            - by: The variables (column names) to group samples.
                Multiple columns should be separated by `,`.
            - method (choice): The method to calculate diversity.
                - gini: The Gini coefficient.
                    It measures the inequality among values of a frequency
                    distribution (for example levels of income).
                - gini.simp: The Gini-Simpson index.
                    It is the probability of interspecific encounter, i.e.,
                    probability that two entities represent different types.
                - inv.simp: Inverse Simpson index.
                    It is the effective number of types that is obtained when
                    the weighted arithmetic mean is used to quantify average
                    proportional abundance of types in the dataset of interest.
                - div: true diversity, or the effective number of types.
                    It refers to the number of equally abundant types needed
                    for the average proportional abundance of the types to
                    equal that observed in the dataset of interest where all
                    types may not be equally abundant.
            - devpars (ns): The parameters for the plotting device.
                - width (type=int): The width of the device
                - height (type=int): The height of the device
                - res (type=int): The resolution of the device
            - cases (type=json): If you have multiple cases, you can use this
                argument to specify them. The keys will be the names of the
                cases. The values will be passed to the corresponding arguments
                above. If any of these arguments are not specified, the values
                in `envs.sample_diversity` will be used. If NO cases are
                specified, the default case will be added, with the name
                `DEFAULT`.

    Requires:
        r-immunarch:
            - check: {{proc.lang}} -e "library(immunarch)"
    """
    input = "immfile:file"
    output = "outdir:dir:{{in.immfile | stem}}.tcrclusters_stats"
    lang = config.lang.rscript
    envs = {
        "cluster_size": {
            "by": "Sample",
            "devpars": {"width": 1000, "height": 900, "res": 100},
            "cases": {},
        },
        "shared_clusters": {
            "numbers_on_heatmap": True,
            "heatmap_meta": [],
            "grouping": None,
            "devpars": {"width": 1000, "height": 1000, "res": 100},
            "cases": {},
        },
        "sample_diversity": {
            "by": None,
            "method": "gini",
            "devpars": {"width": 1000, "height": 1000, "res": 100},
            "cases": {},
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


class CDR3AAPhyschem(Proc):
    """CDR3 AA physicochemical feature analysis

    The idea is to perform a regression between two groups of cells
    (e.g. Treg vs Tconv) at different length of CDR3 AA sequences.
    The regression will be performed for each physicochemical feature of the
    AA (hydrophobicity, volume and isolectric point).

    Reference:
        - https://www.nature.com/articles/ni.3491
        - https://www.nature.com/articles/s41590-022-01129-x
        - Wimley, W. C. & White, S. H. Experimentally determined hydrophobicity
            scale for proteins at membrane - interfaces. Nat. Struct. Biol. 3,
            842-848 (1996).
        - Hdbk of chemistry & physics 72nd edition. (CRC Press, 1991).
        - Zamyatnin, A. A. Protein volume in solution. Prog. Biophys. Mol. Biol.
            24, 107-123 (1972).

    Input:
        immdata: The data loaded by `immunarch::repLoad()`, saved in RDS format
        srtobj: The `Seurat` object, saved in RDS format, used to get the
            metadata for each cell (e.g. cell type)
            It could also be a tab delimited file with `meta.data` of the
            `Seurat` object.
            It has to have a `Sample` column, which is used to match the
            `immdata` object.
            It is optional, if not provided, the metadata from the `immdata`
            object will be used.

    Output:
        outdir: The output directory

    Envs:
        group: The key of group in metadata to define the groups to
            compare. For example, `CellType`, which has cell types annotated
            for each cell in the combined object (immdata + Seurat metadata)
        comparison (type=json): A dict of two groups, with keys as the
            group names and values as the group labels. For example,
            >>> Treg = ["CD4 CTL", "CD4 Naive", "CD4 TCM", "CD4 TEM"]
            >>> Tconv = "Tconv"
        prefix: The prefix of the cell names (rownames) in the metadata.
            The prefix is usually not needed in immdata, as the data is stored
            in the `immdata` object separately for each sample. However, the
            `Seurat` object has a combined `meta.data` for all the samples,
            so the prefix is needed. Usually, the prefix is the sample name.
            For example, `Sample1-AACGTTGAGGCTACGT-1`.
            We need this prefix to add the sample name to the cell names in
            immdata, so that we can match the cells in `immdata` and
            `Seurat` object. Set it to `None` or an empty string if the
            `Seurat` object has the same cell names as `immdata`. You can use
            placeholders to specify the prefix, e.g., `{Sample}_`. In such a
            case, the `Sample` column must exist in the `Seurat` object.
        target: Which group to use as the target group. The target
            group will be labeled as 1, and the other group will be labeled as
            0 in the regression.
        subset: A column, or a list of columns separated by comma,
            in the merged object to subset the cells to perform the regression,
            for each group in the columns.
            If not provided, all the cells will be used.
    """
    input = "immdata:file,srtobj:file"
    output = "outdir:dir:{{in.immdata | stem}}.cdr3aaphyschem"
    lang = config.lang.rscript
    envs = {
        "group": None,
        "comparison": None,
        "prefix": "{Sample}_",
        "target": None,
        "subset": None,
    }
    script = "file://../scripts/tcr/CDR3AAPhyschem.R"
    plugin_opts = {"report": "file://../reports/tcr/CDR3AAPhyschem.svelte"}
