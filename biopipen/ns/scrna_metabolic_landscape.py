"""Metabolic landscape analysis for scRNA-seq data"""
from __future__ import annotations
from pathlib import Path
from typing import Type

from diot import Diot
from datar.tibble import tibble
from pipen_cli_run import Pipeline, process

from ..core.config import config
from ..core.proc import Proc

DEFAULT_OPTS = Diot(
    metafile=None,
    is_seurat=None,
    gmtfile=None,
    grouping=None,
    grouping_prefix="",
    subsetting=None,
    subsetting_prefix=None,
    subsetting_comparison={},
    mutaters=None,
    ncores=config.misc.ncores,
)


class MetabolicPathwayActivity(Proc):
    """Pathway activities for each group

    Requires:
        - name: r-scater
        check: |
            {{proc.lang}} <(echo "library(scater)")
        - name: r-reshape2
        check: |
            {{proc.lang}} <(echo "library(reshape2)")
        - name: r-rcolorbrewer
        check: |
            {{proc.lang}} <(echo "library(RColorBrewer)")
        - name: r-ggplot2
        check: |
            {{proc.lang}} <(echo "library(ggplot2)")
        - name: r-ggprism
        check: |
            {{proc.lang}} <(echo "library(ggprism)")
        - name: r-complexheatmap
        check: |
            {{proc.lang}} <(echo "library(ComplexHeatmap)")
        - name: r-parallel
        check: |
            {{proc.lang}} <(echo "library(parallel)")
    """
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
    """Inter-subset metabolic features - Enrichment analysis in details

    Requires:
        - name: r-parallel
        check: |
            {{proc.lang}} <(echo "library(parallel)")
        - name: r-fgsea
        check: |
            {{proc.lang}} <(echo "library(fgsea)")
    """
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
    script = "file://../scripts/scrna_metabolic_landscape/MetabolicFeatures.R"
    plugin_opts = {
        "report": (
            "file://../reports/"
            "scrna_metabolic_landscape/MetabolicFeatures.svelte"
        )
    }


class MetabolicFeaturesIntraSubset(Proc):
    """Intra-subset metabolic features - Enrichment analysis in details

    Requires:
        - name: r-parallel
        check: |
            {{proc.lang}} <(echo "library(parallel)")
        - name: r-scater
        check: |
            {{proc.lang}} <(echo "library(scater)")
        - name: r-fgsea
        check: |
            {{proc.lang}} <(echo "library(fgsea)")
    """
    input = "sobjfile:file"
    output = "outdir:dir:{{in.sobjfile | stem}}.intra-subset-pathwayfeatures"
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
        "MetabolicFeaturesIntraSubsets.R"
    )
    plugin_opts = {
        "report": (
            "file://../reports/scrna_metabolic_landscape/"
            "MetabolicFeaturesIntraSubsets.svelte"
        )
    }


class MetabolicPathwayHeterogeneity(Proc):
    """Pathway heterogeneity

    Requires:
        - name: r-gtools
        check: |
            {{proc.lang}} <(echo "library(gtools)")
        - name: r-ggplot2
        check: |
            {{proc.lang}} <(echo "library(ggplot2)")
        - name: r-ggprism
        check: |
            {{proc.lang}} <(echo "library(ggprism)")
        - name: r-parallel
        check: |
            {{proc.lang}} <(echo "library(parallel)")
        - name: r-dplyr
        check: |
            {{proc.lang}} <(echo "library(dplyr)")
        - name: r-tibble
        check: |
            {{proc.lang}} <(echo "library(tibble)")
        - name: r-enrichr
        check: |
            {{proc.lang}} <(echo "library(enrichR)")
        - name: r-data.table
        check: |
            {{proc.lang}} <(echo "library(data.table)")
        - name: r-fgsea
        check: |
            {{proc.lang}} <(echo "library(fgsea)")
    """
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


class ScrnaMetabolicLandscape(Pipeline):
    """Metabolic landscape analysis for scRNA-seq data

    An abstract from
    https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape

    See docs here for more details
    https://pwwang.github.io/biopipen/pipelines/scrna_metabolic_landscape


    Reference:
        Xiao, Zhengtao, Ziwei Dai, and Jason W. Locasale.
        "Metabolic landscape of the tumor microenvironment at
        single cell resolution." Nature communications 10.1 (2019): 1-12.

    Input files:
    - metafile: Either a metafile or an rds file of a Seurat object.
        If it is a metafile, it should have two columns: `Sample` and
        `RNADir`. `Sample` should be the first column with unique
        identifiers for the samples and `RNADir` indicates where the
        barcodes, genes, expression matrices are. The data will be loaded
        and an unsupervised clustering will be done.
        Currently only 10X data is supported.
        If it is an rds file, the seurat object will be used directly
    - is_seurat: Whether the input `metafile` is a seurat object.
        If `metafile` is specified directly, this option will be ignored
        and will be inferred from the file extension. If `metafile` is
        not specified, meaning `<pipeline>.procs.MetabolicInput` is dependent
        on other processes, this option will be used to determine whether
        the input is a seurat object or not.
    - gmtfile: The GMT file with the metabolic pathways. The gene names should
        match the gene names in the gene list in RNADir or the Seurat object

    Global options:
    - grouping: defines the basic groups to investigate the metabolic activity
        Typically the clusters.
    - grouping_prefix: Working as a prefix to group names
        For example, if we have `grouping_prefix = "cluster"` and
        we have `1` and `2` in the `grouping` column, the groups
        will be named as `cluster_1` and `cluster_2`
    - subsetting: How do we subset the data. Another column in the metadata
        to do comparisons.
    - subsetting_prefix: Working as a prefix to subset names
        For example, if we have `subsetting_prefix = "timepoint"` and
        we have `pre` and `post` in the `subsetting` column, the subsets
        will be named as `timepoint_pre` and `timepoint_post`
    - subsetting_comparison: What kind of comparisons are we doing to compare
        cells from different subsets.
        It should be dict with keys as the names of the comparisons and
        values as the 2 comparison groups from the `subsetting` column.
        For example, if we have `pre` and `post` in the `subsetting` column,
        we could have `subsetting_comparison = {"pre_vs_post": ["post", "pre"]}`
        The second group will be the control group in the comparison.
        If we also have `1`, `2` and `3` in the `grouping` column, by default,
        the comparisons are done within each subset for each group. For example,
        For group `1`, groups `2` and `3` will be used as control, and for
        group `2`, groups `1` and `3` will be used as control, and for group
        `3`, groups `1` and `2` will be used as control. It is similar to
        `Seurat::FindMarkers` procedure. With this option, the comparisons
        are also done to compare cells from different subsets within each group.
        With the example above, we will have `pre_vs_post` comparisons within
        each group.
    - mutaters: Add new columns to the metadata for grouping/subsetting.
        They are passed to `sobj@meta.data |> mutate(...)`. For example,
        `{"timepoint": "if_else(treatment == 'control', 'pre', 'post')"}`
        will add a new column `timepoint` to the metadata with values of
        `pre` and `post` based on the `treatment` column.
    - ncores: Number of cores to use for parallelization for each process
    """

    defaults = config.pipeline.scrna_metabolic_landscape

    @process(start=True)
    def build_input(self) -> Type[Proc]:
        """Build MetabolicInputs process"""
        from .misc import File2Proc

        if self.options.metafile:
            suffix = Path(self.options.metafile).suffix
            self.options.is_seurat = suffix in (".rds", ".RDS")

        class MetabolicInput(File2Proc):
            """Input for the metabolic pathway analysis pipeline for
            scRNA-seq data

            Input:
                metafile: A metafile indicating the metadata or the rds file
                    with seruat object

            Output:
                metafile: Soft link to `in.metafile`
            """

            if self.options.metafile:
                input_data = [self.options.metafile]

        return MetabolicInput

    @process
    def build_preparing(self, input_proc: Type[Proc]) -> Type[Proc]:
        """Build SeuratPreparing process"""
        from .scrna import SeuratPreparing

        class SeuratPreparing(SeuratPreparing):
            requires = input_proc

        return SeuratPreparing

    @process
    def build_clustering(self, preparing_proc: Type[Proc]) -> Type[Proc]:
        """Build SeuratClustering process"""
        from .scrna import SeuratClustering

        class SeuratClustering(SeuratClustering):
            requires = preparing_proc

        return SeuratClustering

    @process
    def build_mutater(self, clustering_proc: Type[Proc]) -> Type[Proc]:
        """Build SeuratMetadataMutater process"""
        from .scrna import SeuratMetadataMutater

        class SeuratMetadataMutater(SeuratMetadataMutater):
            requires = clustering_proc
            input_data = lambda ch: tibble(
                srtobj=ch.iloc[:, 0],
                metafile=[None],
                mutaters=[self.options.mutaters],
            )

        return SeuratMetadataMutater

    @process
    def build_expr_impute(self, subset_proc: Type[Proc]) -> Type[Proc]:
        """Build MetabolicExprImpute process"""
        from .scrna import ExprImpute

        class MetabolicExprImpute(ExprImpute):
            requires = subset_proc

        return MetabolicExprImpute

    @process(end=True)
    def build_pathway_activity(self, expr_imp_proc: Type[Proc]) -> Type[Proc]:
        """Build MetabolicPathwayActivity process"""
        return Proc.from_proc(
            MetabolicPathwayActivity,
            "MetabolicPathwayActivity",
            requires=expr_imp_proc,
            order=-1,
            envs={
                "ncores": self.options.ncores,
                "gmtfile": self.options.gmtfile,
                "grouping": self.options.grouping,
                "grouping_prefix": self.options.grouping_prefix,
                "subsetting": self.options.subsetting,
                "subsetting_prefix": self.options.subsetting_prefix,
            },
        )

    @process(end=True)
    def build_pathway_heterogeneity(self, norm_proc: Type[Proc]) -> Type[Proc]:
        """Build MetabolicPathwayHeterogeneity process"""
        return Proc.from_proc(
            MetabolicPathwayHeterogeneity,
            "MetabolicPathwayHeterogeneity",
            requires=norm_proc,
            envs={
                "ncores": self.options.ncores,
                "gmtfile": self.options.gmtfile,
                "grouping": self.options.grouping,
                "grouping_prefix": self.options.grouping_prefix,
                "subsetting": self.options.subsetting,
                "subsetting_prefix": self.options.subsetting_prefix,
            },
        )

    @process(end=True)
    def build_features(self, norm_proc: Type[Proc]) -> Type[Proc]:
        """Build MetabolicFeatures process"""
        return Proc.from_proc(
            MetabolicFeatures,
            "MetabolicFeatures",
            requires=norm_proc,
            envs={
                "ncores": self.options.ncores,
                "gmtfile": self.options.gmtfile,
                "grouping": self.options.grouping,
                "grouping_prefix": self.options.grouping_prefix,
                "subsetting": self.options.subsetting,
                "subsetting_prefix": self.options.subsetting_prefix,
            },
        )

    @process(end=True)
    def build_features_intra_subset(self, norm_proc: Type[Proc]) -> Type[Proc]:
        """Build MetabolicFeaturesIntraSubset process"""
        if self.options.subsetting_comparison and not self.options.subsetting:
            raise ValueError(
                "Cannot use `subsetting_comparison` without `subsetting`."
            )

        return Proc.from_proc(
            MetabolicFeaturesIntraSubset,
            "MetabolicFeaturesIntraSubset",
            requires=norm_proc,
            envs={
                "ncores": self.options.ncores,
                "gmtfile": self.options.gmtfile,
                "grouping": self.options.grouping,
                "grouping_prefix": self.options.grouping_prefix,
                "subsetting": self.options.subsetting,
                "subsetting_prefix": self.options.subsetting_prefix,
                "subsetting_comparison": self.options.subsetting_comparison,
            },
        )

    def build(self) -> None:
        """Build processes for metabolic landscape analysis pipeline"""
        self.options = DEFAULT_OPTS | self.options

        if not self.options.gmtfile:
            raise ValueError("`gmtfile` with metabolic pathways is required.")

        MetabolicInput = self.build_input()

        if self.options.is_seurat:
            # Use the rds file
            SeuratClustering = MetabolicInput
        else:
            # Do clustering
            SeuratPreparing = self.build_preparing(MetabolicInput)
            SeuratClustering = self.build_clustering(SeuratPreparing)

        if self.options.mutaters:
            SeuratMetadataMutater = self.build_mutater(SeuratClustering)
        else:
            # No mutaters, just use the SeuratClustering
            SeuratMetadataMutater = SeuratClustering

        # Do imputation and normalization for all the data together
        MetabolicExprImpute = self.build_expr_impute(SeuratMetadataMutater)

        self.build_pathway_activity(MetabolicExprImpute)
        self.build_pathway_heterogeneity(MetabolicExprImpute)
        self.build_features(MetabolicExprImpute)
        if self.options.subsetting_comparison:
            self.build_features_intra_subset(MetabolicExprImpute)
