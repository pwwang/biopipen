"""Metabolic landscape analysis for scRNA-seq data"""
from __future__ import annotations
from pathlib import Path
from typing import Type

from diot import Diot
from datar.tibble import tibble
from pipen.utils import mark
from pipen_args import ProcGroup
from pipen_annotate import annotate

from ..core.config import config
from ..core.proc import Proc


class MetabolicPathwayActivity(Proc):
    """This process calculates the pathway activities in different groups and subsets.

    The cells are first grouped by subsets and then the metabolic activities are
    examined for each groups in different subsets.

    For each subset, a heatmap and a violin plot will be generated.
    The heatmap shows the pathway activities for each group and each metabolic pathway

    ![MetabolicPathwayActivity_heatmap](https://pwwang.github.io/immunopipe/latest/processes/images/MetabolicPathwayActivity_heatmap.png){: width="80%"}

    The violin plot shows the distribution of the pathway activities for each group

    ![MetabolicPathwayActivity_violin](https://pwwang.github.io/immunopipe/latest/processes/images/MetabolicPathwayActivity_violin.png){: width="45%"}

    Envs:
        ntimes (type=int): Number of times to do the permutation
        ncores (type=int;pgarg): Number of cores to use for parallelization
            Defaults to `ScrnaMetabolicLandscape.ncores`
        heatmap_devpars (ns): Device parameters for the heatmap
            - width (type=int): Width of the heatmap
            - height (type=int): Height of the heatmap
            - res (type=int): Resolution of the heatmap
        violin_devpars (ns): Device parameters for the violin plot
            - width (type=int): Width of the violin plot
            - height (type=int): Height of the violin plot
            - res (type=int): Resolution of the violin plot
        gmtfile (pgarg): The GMT file with the metabolic pathways.
            Defaults to `ScrnaMetabolicLandscape.gmtfile`
        grouping (type=auto;pgarg;readonly): Defines the basic groups to
            investigate the metabolic activity, typically the clusters.
            Defaults to `ScrnaMetabolicLandscape.grouping`
        grouping_prefix (type=auto;pgarg;readonly): Working as a prefix to group
            names. For example, if we have `grouping_prefix = "cluster"` and
            we have `1` and `2` in the `grouping` column, the groups
            will be named as `cluster_1` and `cluster_2`.
            Defaults to `ScrnaMetabolicLandscape.grouping_prefix`
        subsetting (type=auto;pgarg;readonly): How do we subset the data. Other
            columns in the metadata to do comparisons. For example,
            `"TimePoint"` or `["TimePoint", "Response"]`.
            Defaults to `ScrnaMetabolicLandscape.subsetting`
        subsetting_prefix (type=auto;pgarg;readonly): Working as a prefix to
            subset names.
            For example, if we have `subsetting_prefix = "timepoint"` and
            we have `pre` and `post` in the `subsetting` column, the subsets
            will be named as `timepoint_pre` and `timepoint_post`.
            If `subsetting` is a list, then this should also be a
            same-length list. If a single string is given, it will be
            repeated to a list with the same length as `subsetting`.
            Defaults to `ScrnaMetabolicLandscape.subsetting_prefix`

    Requires:
        r-scater:
            - check: {{proc.lang}} <(echo "library(scater)")
        r-reshape2:
            - check: {{proc.lang}} <(echo "library(reshape2)")
        r-rcolorbrewer:
            - check: {{proc.lang}} <(echo "library(RColorBrewer)")
        r-ggplot2:
            - check: {{proc.lang}} <(echo "library(ggplot2)")
        r-ggprism:
            - check: {{proc.lang}} <(echo "library(ggprism)")
        r-complexheatmap:
            - check: {{proc.lang}} <(echo "library(ComplexHeatmap)")
        r-parallel:
            - check: {{proc.lang}} <(echo "library(parallel)")
    """  # noqa: E501
    input = "sobjfile:file"
    output = "outdir:dir:{{in.sobjfile | stem}}.pathwayactivity"
    envs = {
        "ntimes": 5000,
        "ncores": config.misc.ncores,
        "heatmap_devpars": {},
        "violin_devpars": {},
        "gmtfile": None,
        "grouping": None,
        "grouping_prefix": "",
        "subsetting": None,
        "subsetting_prefix": "",
    }
    lang = config.lang.rscript
    script = (
        "file://../scripts/"
        "scrna_metabolic_landscape/MetabolicPathwayActivity.R"
    )
    plugin_opts = {
        "report": (
            "file://../reports/"
            "scrna_metabolic_landscape/MetabolicPathwayActivity.svelte"
        )
    }


class MetabolicFeatures(Proc):
    """This process performs enrichment analysis for the metabolic pathways
    for each group in each subset.

    The enrichment analysis is done with [`fgsea`](https://bioconductor.org/packages/release/bioc/html/fgsea.html)
    package or the [`GSEA_R`](https://github.com/GSEA-MSigDB/GSEA_R) package.

    Envs:
        ncores (type=int;pgarg): Number of cores to use for parallelization.
            Defaults to `ScrnaMetabolicLandscape.ncores`
        fgsea (flag): Whether to do fast gsea analysis using `fgsea` package.
            If `False`, the `GSEA_R` package will be used.
        prerank_method (choice): Method to use for gene preranking.
            Signal to noise: the larger the differences of the means
            (scaled by the standard deviations); that is, the more distinct
            the gene expression is in each phenotype and the more the gene
            acts as a “class marker.”.
            Absolute signal to noise: the absolute value of the signal to
            noise.
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
        top (type=int): N top of enriched pathways to show
        gmtfile (pgarg): The GMT file with the metabolic pathways.
            Defaults to `ScrnaMetabolicLandscape.gmtfile`
        grouping (type=auto;pgarg;readonly): Defines the basic groups to
            investigate the metabolic activity.
            Defaults to `ScrnaMetabolicLandscape.grouping`
        grouping_prefix (type=auto;pgarg;readonly): Working as a prefix to
            group names.
            Defaults to `ScrnaMetabolicLandscape.grouping_prefix`
        subsetting (type=auto;pgarg;readonly): How do we subset the data.
            Another column(s) in the metadata.
            Defaults to `ScrnaMetabolicLandscape.subsetting`
        subsetting_prefix (type=auto;pgarg;readonly): Working as a prefix to
            subset names.
            Defaults to `ScrnaMetabolicLandscape.subsetting_prefix`

    Requires:
        r-parallel:
            - check: {{proc.lang}} <(echo "library(parallel)")
        r-fgsea:
            - check: {{proc.lang}} <(echo "library(fgsea)")
    """  # noqa: E501
    input = "sobjfile:file"
    output = "outdir:dir:{{in.sobjfile | stem}}.pathwayfeatures"
    lang = config.lang.rscript
    envs = {
        "ncores": config.misc.ncores,
        "fgsea": True,
        "prerank_method": "signal_to_noise",
        "top": 10,
        "gmtfile": None,
        "grouping": None,
        "grouping_prefix": "",
        "subsetting": None,
        "subsetting_prefix": "",
    }
    script = (
        "file://../scripts/scrna_metabolic_landscape/MetabolicFeatures.R"
    )
    plugin_opts = {
        "report": (
            "file://../reports/"
            "scrna_metabolic_landscape/MetabolicFeatures.svelte"
        )
    }


class MetabolicFeaturesIntraSubset(Proc):
    """Intra-subset metabolic features - Enrichment analysis in details

    Similar to the [`MetabolicFeatures`](!!#biopipennsscrna_metabolic_landscapemetabolicfeatures)
    process, this process performs enrichment analysis for the metabolic pathways for
    each subset in each group, instead of each group in each subset.

    Envs:
        ncores (type=int; pgarg): Number of cores to use for parallelization
            Defaults to `ScrnaMetabolicLandscape.ncores`
        fgsea (flag): Whether to do fast gsea analysis
        prerank_method (choice): Method to use for gene preranking
            Signal to noise: the larger the differences of the means
            (scaled by the standard deviations); that is, the more distinct
            the gene expression is in each phenotype and the more the gene
            acts as a “class marker.”.
            Absolute signal to noise: the absolute value of the signal to
            noise.
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
        top (type=int): N top of enriched pathways to show
        gmtfile (pgarg): The GMT file with the metabolic pathways.
            Defaults to `ScrnaMetabolicLandscape.gmtfile`
        grouping (type=auto;pgarg;readonly): Defines the basic groups to
            investigate the metabolic activity.
            Defaults to `ScrnaMetabolicLandscape.grouping`
        grouping_prefix (type=auto;pgarg;readonly): Working as a prefix to group
            names.
            Defaults to `ScrnaMetabolicLandscape.grouping_prefix`
        subsetting (type=auto;pgarg;readonly): How do we subset the data.
            Another column(s) in the metadata.
            Defaults to `ScrnaMetabolicLandscape.subsetting`
        subsetting_prefix (type=auto;pgarg;readonly): Working as a prefix to
            subset names.
            Defaults to `ScrnaMetabolicLandscape.subsetting_prefix`
        subsetting_comparison (type=json;pgarg;readonly): How do we compare the
            subsets.
            Defaults to `ScrnaMetabolicLandscape.subsetting_comparison`

    Requires:
        r-parallel:
            - check: {{proc.lang}} <(echo "library(parallel)")
        r-scater:
            - check: {{proc.lang}} <(echo "library(scater)")
        r-fgsea:
            - check: {{proc.lang}} <(echo "library(fgsea)")
    """  # noqa: E501
    input = "sobjfile:file"
    output = (
        "outdir:dir:{{in.sobjfile | stem}}.intra-subset-pathwayfeatures"
    )
    lang = config.lang.rscript
    envs = {
        "ncores": config.misc.ncores,
        "gmtfile": None,
        "fgsea": True,
        "prerank_method": "signal_to_noise",
        "top": 10,
        "grouping": None,
        "grouping_prefix": "",
        "subsetting": None,
        "subsetting_prefix": "",
        "subsetting_comparison": {},
    }
    script = (
        "file://../scripts/scrna_metabolic_landscape/"
        "MetabolicFeaturesIntraSubset.R"
    )
    plugin_opts = {
        "report": (
            "file://../reports/scrna_metabolic_landscape/"
            "MetabolicFeaturesIntraSubset.svelte"
        )
    }


class MetabolicPathwayHeterogeneity(Proc):
    """Calculate Metabolic Pathway heterogeneity.

    For each subset, the normalized enrichment score (NES) of each metabolic pathway
    is calculated for each group.
    The NES is calculated by comparing the enrichment score of the subset to the
    enrichment scores of the same subset in the permutations.
    The p-value is calculated by comparing the NES to the NESs of the same subset
    in the permutations.
    The heterogeneity can be reflected by the NES values and the p-values in
    different groups for the metabolic pathways.

    ![MetabolicPathwayHeterogeneity](https://pwwang.github.io/immunopipe/latest/processes/images/MetabolicPathwayHeterogeneity.png)


    Envs:
        gmtfile (pgarg): The GMT file with the metabolic pathways.
            Defaults to `ScrnaMetabolicLandscape.gmtfile`
        select_pcs (type=float): Select the PCs to use for the analysis.
        pathway_pval_cutoff (type=float): The p-value cutoff to select
            the enriched pathways
        ncores (type=int;pgarg): Number of cores to use for parallelization
            Defaults to `ScrnaMetabolicLandscape.ncores`
        bubble_devpars (ns): The devpars for the bubble plot
            - width (type=int): The width of the plot
            - height (type=int): The height of the plot
            - res (type=int): The resolution of the plot
        grouping (type=auto;pgarg;readonly): Defines the basic groups to
            investigate the metabolic activity.
            Defaults to `ScrnaMetabolicLandscape.grouping`
        grouping_prefix (type=auto;pgarg;readonly): Working as a prefix to group
            names.
            Defaults to `ScrnaMetabolicLandscape.grouping_prefix`
        subsetting (type=auto;pgarg;readonly): How do we subset the data.
            Another column(s) in the metadata.
            Defaults to `ScrnaMetabolicLandscape.subsetting`
        subsetting_prefix (type=auto;pgarg;readonly): Working as a prefix to
            subset names.
            Defaults to `ScrnaMetabolicLandscape.subsetting_prefix`

    Requires:
        r-gtools:
            - check: {{proc.lang}} <(echo "library(gtools)")
        r-ggplot2:
            - check: {{proc.lang}} <(echo "library(ggplot2)")
        r-ggprism:
            - check: {{proc.lang}} <(echo "library(ggprism)")
        r-parallel:
            - check: {{proc.lang}} <(echo "library(parallel)")
        r-dplyr:
            - check: {{proc.lang}} <(echo "library(dplyr)")
        r-tibble:
            - check: {{proc.lang}} <(echo "library(tibble)")
        r-enrichr:
            - check: {{proc.lang}} <(echo "library(enrichR)")
        r-data.table:
            - check: {{proc.lang}} <(echo "library(data.table)")
        r-fgsea:
            - check: {{proc.lang}} <(echo "library(fgsea)")
    """  # noqa: E501
    input = "sobjfile:file"
    output = "outdir:dir:{{in.sobjfile | stem}}.pathwayhetero"
    lang = config.lang.rscript
    envs = {
        "gmtfile": None,
        "select_pcs": 0.8,
        "pathway_pval_cutoff": 0.01,
        "ncores": config.misc.ncores,
        "bubble_devpars": {},
        "grouping": None,
        "grouping_prefix": "",
        "subsetting": None,
        "subsetting_prefix": "",
    }
    script = (
        "file://../scripts/scrna_metabolic_landscape/"
        "MetabolicPathwayHeterogeneity.R"
    )
    plugin_opts = {
        "report": (
            "file://../reports/scrna_metabolic_landscape/"
            "MetabolicPathwayHeterogeneity.svelte"
        )
    }


class ScrnaMetabolicLandscape(ProcGroup):
    """Metabolic landscape analysis for scRNA-seq data

    An abstract from
    <https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape>

    See docs here for more details
    <https://pwwang.github.io/biopipen/pipelines/scrna_metabolic_landscape>


    Reference:
        Xiao, Zhengtao, Ziwei Dai, and Jason W. Locasale.
        "Metabolic landscape of the tumor microenvironment at
        single cell resolution." Nature communications 10.1 (2019): 1-12.

    Args:
        metafile: Either a metafile or an rds file of a Seurat object.
            If it is a metafile, it should have two columns: `Sample` and
            `RNAData`. `Sample` should be the first column with unique
            identifiers for the samples and `RNAData` indicates where the
            barcodes, genes, expression matrices are. The data will be loaded
            and an unsupervised clustering will be done.
            Currently only 10X data is supported.
            If it is an rds file, the seurat object will be used directly
        is_seurat (flag): Whether the input `metafile` is a seurat object.
            If `metafile` is specified directly, this option will be ignored
            and will be inferred from the file extension. If `metafile` is
            not specified, meaning `<pipeline>.procs.MetabolicInput` is
            dependent on other processes, this option will be used to determine
            whether the input is a seurat object or not.
        noimpute (flag): Whether to do imputation for the dropouts.
            If False, the values will be left as is.
        gmtfile: The GMT file with the metabolic pathways. The gene names should
            match the gene names in the gene list in RNAData or
            the Seurat object.
            You can also provide a URL to the GMT file.
            For example, from
            <https://download.baderlab.org/EM_Genesets/current_release/Human/symbol/>.
        grouping: defines the basic groups to investigate the metabolic activity
            Typically the clusters.
        grouping_prefix: Working as a prefix to group names
            For example, if we have `grouping_prefix = "cluster"` and
            we have `1` and `2` in the `grouping` column, the groups
            will be named as `cluster_1` and `cluster_2`
        subsetting (type=auto): How do we subset the data. Other columns in the
            metadata to do comparisons. For example, `"TimePoint"` or
            `["TimePoint", "Response"]`
        subsetting_prefix (type=auto): Working as a prefix to subset names
            For example, if we have `subsetting_prefix = "timepoint"` and
            we have `pre` and `post` in the `subsetting` column, the subsets
            will be named as `timepoint_pre` and `timepoint_post`
            If `subsetting` is a list, then this should also be a same-length
            list. If a single string is given, it will be repeated to a list
            with the same length as `subsetting`
        subsetting_comparison (type=json): What kind of comparisons are we
            doing to compare cells from different subsets.
            It should be dict with keys as the names of the comparisons and
            values as the 2 comparison groups from the `subsetting` column.
            For example, if we have `pre` and `post` in the `subsetting` column,
            we could have
            `subsetting_comparison = {"pre_vs_post": ["post", "pre"]}`
            The second group will be the control group in the comparison.
            If we also have `1`, `2` and `3` in the `grouping` column,
            by default, the comparisons are done within each subset for
            each group. For example, for group `1`, groups `2` and `3`
            will be used as control, and for group `2`, groups `1` and `3`
            will be used as control, and for group `3`, groups `1` and `2`
            will be used as control. It is similar to `Seurat::FindMarkers`
            procedure. With this option, the comparisons are also done to
            compare cells from different subsets within each group. With the
            example above, we will have `pre_vs_post` comparisons within
            each group.
            If `subsetting` is a list, this must be a list of dicts with the
            same length.
        mutaters (type=json): Add new columns to the metadata for
            grouping/subsetting.
            They are passed to `sobj@meta.data |> mutate(...)`. For example,
            `{"timepoint": "if_else(treatment == 'control', 'pre', 'post')"}`
            will add a new column `timepoint` to the metadata with values of
            `pre` and `post` based on the `treatment` column.
        ncores (type=int): Number of cores to use for parallelization for
            each process
    """
    DEFAULTS = Diot(
        metafile=None,
        is_seurat=None,
        gmtfile=None,
        grouping=None,
        grouping_prefix="",
        subsetting=None,
        subsetting_prefix=None,
        subsetting_comparison={},
        mutaters=None,
        noimpute=False,
        ncores=config.misc.ncores,
    )

    def post_init(self):
        """Load runtime processes"""
        if self.opts.metafile:
            suffix = Path(self.opts.metafile).suffix
            self.opts.is_seurat = suffix in (".rds", ".RDS")

        # Make sure the grouping is a list
        if self.opts.subsetting and not isinstance(self.opts.subsetting, list):
            self.opts.subsetting = [self.opts.subsetting]

        # Make sure the grouping is a list with the same length as subsetting
        if (
            self.opts.subsetting
            and not isinstance(self.opts.subsetting_prefix, list)
        ):
            self.opts.subsetting_prefix = [
                self.opts.subsetting_prefix
            ] * len(self.opts.subsetting)

        # Make sure the lengths of subsetting and subsetting_comparison the same
        if self.opts.subsetting:
            if len(self.opts.subsetting) == 1 and isinstance(
                self.opts.subsetting_comparison, dict
            ):
                self.opts.subsetting_comparison = [
                    self.opts.subsetting_comparison
                ]

            if len(self.opts.subsetting) > 1 and not isinstance(
                self.opts.subsetting_comparison, list
            ):
                raise ValueError(
                    "The length of `subsetting` is larger than 1, "
                    "but `subsetting_comparison` is not a list of dicts."
                )

            if len(self.opts.subsetting) != len(
                self.opts.subsetting_comparison
            ):
                raise ValueError(
                    "The length of `subsetting` and `subsetting_comparison` "
                    "are not the same"
                )

    @ProcGroup.add_proc
    def p_input(self) -> Type[Proc]:
        """Build MetabolicInputs process"""
        from .misc import File2Proc

        @mark(board_config_hidden=True)
        class MetabolicInput(File2Proc):
            """This process takes Seurat object as input and pass it to the next
            processes in the `ScrnaMetabolicLandscape` group.

            There is no configuration for this process.
            """

            if self.opts.metafile:
                input_data = [self.opts.metafile]

        return MetabolicInput

    @ProcGroup.add_proc
    def p_preparing(self) -> Type[Proc]:
        """Build SeuratPreparing process"""
        if self.opts.is_seurat:
            return None

        from .scrna import SeuratPreparing

        class MetabolicSeuratPreparing(SeuratPreparing):
            requires = self.p_input

        return MetabolicSeuratPreparing

    @ProcGroup.add_proc
    def p_clustering(self) -> Type[Proc]:
        """Build SeuratClustering process"""
        if self.opts.is_seurat:
            return self.p_input

        from .scrna import SeuratClustering

        class MetabolicSeuratClustering(SeuratClustering):
            requires = self.p_preparing

        return MetabolicSeuratClustering

    @ProcGroup.add_proc
    def p_mutater(self) -> Type[Proc]:
        """Build SeuratMetadataMutater process"""
        if not self.opts.mutaters:
            return self.p_clustering

        from .scrna import SeuratMetadataMutater

        class MetabolicSeuratMetadataMutater(SeuratMetadataMutater):
            requires = self.p_clustering
            input_data = lambda ch: tibble(
                srtobj=ch.iloc[:, 0],
                metafile=[None],
            )
            envs = {"mutaters": self.opts.mutaters}

        return MetabolicSeuratMetadataMutater

    @ProcGroup.add_proc
    def p_expr_impute(self) -> Type[Proc]:
        """Build  process"""
        if self.opts.noimpute:
            return self.p_mutater

        from .scrna import ExprImputation

        @annotate.format_doc(indent=3)
        class MetabolicExprImputation(ExprImputation):
            """{{Summary}}

            You can turn off the imputation by setting the `noimpute` option
            of the process group to `True`.
            """
            requires = self.p_mutater

        return MetabolicExprImputation

    @ProcGroup.add_proc
    def p_pathway_activity(self) -> Type[Proc]:
        """Build MetabolicPathwayActivity process"""
        return Proc.from_proc(
            MetabolicPathwayActivity,
            "MetabolicPathwayActivity",
            requires=self.p_expr_impute,
            order=-1,
            envs={
                "ncores": self.opts.ncores,
                "gmtfile": self.opts.gmtfile,
                "grouping": self.opts.grouping,
                "grouping_prefix": self.opts.grouping_prefix,
                "subsetting": self.opts.subsetting,
                "subsetting_prefix": self.opts.subsetting_prefix,
            },
        )

    @ProcGroup.add_proc
    def p_pathway_heterogeneity(self) -> Type[Proc]:
        """Build MetabolicPathwayHeterogeneity process"""
        return Proc.from_proc(
            MetabolicPathwayHeterogeneity,
            "MetabolicPathwayHeterogeneity",
            requires=self.p_expr_impute,
            envs={
                "ncores": self.opts.ncores,
                "gmtfile": self.opts.gmtfile,
                "grouping": self.opts.grouping,
                "grouping_prefix": self.opts.grouping_prefix,
                "subsetting": self.opts.subsetting,
                "subsetting_prefix": self.opts.subsetting_prefix,
            },
        )

    @ProcGroup.add_proc
    def p_features(self) -> Type[Proc]:
        """Build MetabolicFeatures process"""
        return Proc.from_proc(
            MetabolicFeatures,
            "MetabolicFeatures",
            requires=self.p_expr_impute,
            envs={
                "ncores": self.opts.ncores,
                "gmtfile": self.opts.gmtfile,
                "grouping": self.opts.grouping,
                "grouping_prefix": self.opts.grouping_prefix,
                "subsetting": self.opts.subsetting,
                "subsetting_prefix": self.opts.subsetting_prefix,
            },
        )

    @ProcGroup.add_proc
    def p_features_intra_subset(self) -> Type[Proc]:
        """Build MetabolicFeaturesIntraSubset process"""
        if self.opts.subsetting_comparison and not self.opts.subsetting:
            raise ValueError(
                "Cannot use `subsetting_comparison` without `subsetting`."
            )

        return Proc.from_proc(
            MetabolicFeaturesIntraSubset,
            "MetabolicFeaturesIntraSubset",
            requires=self.p_expr_impute,
            envs={
                "ncores": self.opts.ncores,
                "gmtfile": self.opts.gmtfile,
                "grouping": self.opts.grouping,
                "grouping_prefix": self.opts.grouping_prefix,
                "subsetting": self.opts.subsetting,
                "subsetting_prefix": self.opts.subsetting_prefix,
                "subsetting_comparison": self.opts.subsetting_comparison,
            },
        )


if __name__ == "__main__":
    ScrnaMetabolicLandscape().as_pipen().run()
