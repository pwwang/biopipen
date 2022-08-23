"""Metabolic landscape analysis for scRNA-seq data

An abstract from https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape

If you have clustering done somewhere else, you could use `replace_clustering()`

See docs here

https://pwwang.github.io/biopipen/pipelines/scrna_metabolic


Reference:
    [1] Xiao, Zhengtao, Ziwei Dai, and Jason W. Locasale.
    "Metabolic landscape of the tumor microenvironment at
    single cell resolution." Nature communications 10.1 (2019): 1-12.

Start Process:
    MetabolicInputs
"""
from typing import Any, Mapping
from pathlib import Path

from diot import Diot
from datar.tibble import tibble
from pipen import Pipen
from pipen.channel import expand_dir
from ..core.config import config
from ..core.proc import Proc

DEFAULT = {
    "clustered": False,
    "intra-subset": True,
}


def _as_config(conf):
    from pipen_filters.filters import FILTERS
    return FILTERS["config"](conf, loader="toml")


def build_processes(options: Mapping[str, Any] = None) -> Proc:
    """Build processes for metabolic landscape analysis pipeline"""

    from .scrna import (
        ExprImpute,
        SeuratPreparing,
        SeuratSubset,
        SeuratClustering,
        SeuratMetadataMutater,
    )
    options = (
        Diot(DEFAULT)
        | (
            config
            .get("pipeline", {})
            .get("scrna_metabolic", {})
        )
        | (options or {})
    )

    class MetabolicInputs(Proc):
        """Input for the metabolic pathway analysis pipeline for
        scRNA-seq data

        Input:
            metafile: A metafile indicating the metadata or the rds file with
                seruat object if `config.pipeline.scrna_metabolic.clustered` is
                `True`.
                For a meta file, Two columns are required: `Sample` and `RNADir`
                `Sample` should be the first column with unique identifiers
                for the samples and `RNADir` indicates where the expression
                matrices are.
                Currently only 10X data is supported
            gmtfile: The GMT file with the metabolic pathways
            config: The configuration file, string in TOML format or a python
                dictionary as config for the analysis
                (based on `envs.config_fmt`)
                They keys include:
                - name: (optional) Used in reports
                - grouping: How do we group the cells
                  groupby - The column used to group by if it exists
                  mutaters - Add new columns to the metadata to group by
                    They are passed to `sobj@meta.data |> mutate(...)`
                - subsetting: How do we subset the data. The imputation
                  will be done in each subset separately
                  groupby - The column used to subset if it exists
                  alias - The alias of the subset working as a prefix to subset
                    names
                  mutaters - Add new columns to the metadata to subset by
                - design: What kind of comparisons are we doing?
                  It should be the values of subsetting `groupby`s

        Output:
            metafile: Soft link to `in.metafile`
            gmtfile: Soft link to `in.gmtfile`
            configfile: The config file with `in.config` saved

        Requires:
            - name: rtoml
              check: |
                {{proc.lang}} -c "import rtoml"
        """

        input = "metafile:file, gmtfile:file, config:var"
        output = [
            "metafile:file:{{in.metafile | basename}}",
            "gmtfile:file:{{in.gmtfile | basename}}",
            "configfile:file:config.toml",
        ]
        lang = config.lang.python
        script = "file://../scripts/scrna_metabolic/MetabolicInputs.py"


    class MetabolicSeuratPreparing(SeuratPreparing):
        if not options["clustered"]:
            requires = MetabolicInputs


    class MetabolicSeuratClustering(SeuratClustering):
        requires = MetabolicSeuratPreparing


    class MetabolicCellGroups(SeuratMetadataMutater):
        """Group cells for metabolic landscape analysis"""
        if options["clustered"]:
            requires = MetabolicInputs
            input_data = lambda ch: tibble(
                srtobj=ch.metafile,
                metafile=[None],
                mutaters=[
                    _as_config(configfile)["grouping"].get("mutaters", {})
                    for configfile in ch.configfile
                ],
            )
        else:
            requires = MetabolicSeuratClustering, MetabolicInputs
            input_data = lambda ch1, ch2: tibble(
                srtobj=ch1.rdsfile,
                metafile=[None],
                mutaters=[
                    _as_config(configfile)["grouping"].get("mutaters", {})
                    for configfile in ch2.configfile
                ],
            )

    def _get_subset_config(configfiles):
        out = []
        for i, conf in enumerate(configfiles):
            conf = _as_config(conf)
            alias = conf["subsetting"].get("alias", "subset")
            key = f"{i}.{alias}"
            out.append({key: conf["subsetting"]})
        return out

    class MetabolicCellSubsets(SeuratSubset):
        requires = MetabolicCellGroups, MetabolicInputs
        input_data = lambda ch1, ch2: tibble(
            srtobj=ch1.rdsfile,
            subsets=_get_subset_config(ch2.configfile),
        )

    class MetabolicExprImputation(ExprImpute):
        """Impute the dropout values in scRNA-seq data."""
        #                      jobname.case_xxx.RDS
        requires = MetabolicCellSubsets
        input_data = lambda ch: tibble(
            infile=sum(
                (
                    expand_dir(
                        ch.iloc[i:i + 1, :],
                        pattern="*.RDS",
                    ).iloc[:, 0].to_list()
                    for i in range(ch.shape[0])
                ),
                [],
            )
        )

    def _group_imputed_files(impfiles):
        """Use the index in the name to group the imputed files"""
        out = {}
        for impfile in impfiles:
            idx, _ = Path(impfile).stem.split(".", 1)
            idx = int(idx)
            out.setdefault(idx, []).append(impfile)

        return [out[key] for key in sorted(out)]

    class MetabolicPrepareSCE(Proc):
        """Prepare SingleCellExperiment objects

        Requires:
            - name: r-scater
              check: |
                {{proc.lang}} <(echo "library(scater)")
            - name: r-seurat
              check: |
                {{proc.lang}} <(echo "library(Seurat)")
        """

        requires = MetabolicExprImputation, MetabolicInputs
        input = "impfiles:files, gmtfile:file"
        input_data = lambda ch1, ch2: tibble(
            impfiles=_group_imputed_files(ch1.outfile),
            gmtfile=ch2.gmtfile,
        )
        output = (
            "outfile:file:"
            "{{in.impfiles | first | stem | split: '.' | first}}.sce.RDS"
        )
        lang = config.lang.rscript
        envs = {"refexon": config.ref.refexon}
        script = "file://../scripts/scrna_metabolic/MetabolicPrepareSCE.R"


    class MetabolicExprNormalization(Proc):
        """Normalize the expression data using deconvolution

        Requires:
            - name: r-scran
              check: |
                {{proc.lang}} <(echo "library(scran)")
        """

        requires = MetabolicPrepareSCE, MetabolicInputs
        input = "sceobj:file, configfile:file"
        output = "outfile:file:{{in.sceobj | stem0}}.sce.RDS"
        input_data = lambda ch1, ch2: tibble(
            sceobj=ch1.outfile,
            configfile=ch2.configfile,
        )
        envs = {"dropout": 0.75, "refexon": config.ref.refexon}
        lang = config.lang.rscript
        script = (
            "file://../scripts/scrna_metabolic/MetabolicExprNormalization.R"
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

        requires = MetabolicExprNormalization, MetabolicInputs
        input = "sceobj:file, gmtfile:file, configfile:file"
        input_data = lambda ch1, ch2: tibble(
            sceobj=ch1.outfile,
            gmtfile=ch2.gmtfile,
            configfile=ch2.configfile,
        )
        output = "outdir:dir:{{in.sceobj | stem}}.pathwayactivity"
        envs = {
            "ntimes": 5000,
            "ncores": config.misc.ncores,
            "heatmap_devpars": {"res": 100, "width": 800, "height": 400},
            "violin_devpars": {"res": 100, "width": 500, "height": 550},
        }
        order = 1
        lang = config.lang.rscript
        script = "file://../scripts/scrna_metabolic/MetabolicPathwayActivity.R"
        plugin_opts = {
            "report": (
                "file://"
                "../reports/scrna_metabolic/MetabolicPathwayActivity.svelte"
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

        requires = MetabolicExprNormalization, MetabolicInputs
        input = "sceobj:file, gmtfile:file, configfile:file"
        input_data = lambda ch1, ch2: tibble(
            sceobj=ch1.outfile,
            gmtfile=ch2.gmtfile,
            configfile=ch2.configfile,
        )
        output = "outdir:dir:{{in.sceobj | stem}}.pathwayhetero"
        lang = config.lang.rscript
        order = 2
        envs = {
            "select_pcs": 0.8,
            "pathway_pval_cutoff": 0.01,
            "ncores": config.misc.ncores,
            "bubble_devpars": {"res": 100, "width": 1200, "height": 700}
        }
        script = (
            "file://../scripts/scrna_metabolic/MetabolicPathwayHeterogeneity.R"
        )
        plugin_opts = {
            "report": (
                "file://../reports/scrna_metabolic/"
                "MetabolicPathwayHeterogeneity.svelte"
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

        requires = MetabolicExprNormalization, MetabolicInputs
        input = "sceobj:file, gmtfile:file, configfile:file"
        input_data = lambda ch1, ch2: tibble(
            sceobj=ch1.outfile,
            gmtfile=ch2.gmtfile,
            configfile=ch2.configfile,
        )
        output = "outdir:dir:{{in.sceobj | stem}}.pathwayfeatures"
        lang = config.lang.rscript
        order = 3
        envs = {
            "ncores": config.misc.ncores,
            "fgsea": True,
            "prerank_method": "signal_to_noise",
            "top": 10,
        }
        script = "file://../scripts/scrna_metabolic/MetabolicFeatures.R"
        plugin_opts = {
            "report": (
                "file://../reports/scrna_metabolic/MetabolicFeatures.svelte"
            )
        }


    class MetabolicFeaturesIntraSubsets(Proc):
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
        if options["intra-subset"]:
            requires = MetabolicExprNormalization, MetabolicInputs
        input = "sceobjs:files, gmtfile:file, configfile:file"
        input_data = lambda ch1, ch2: tibble(
            sceobjs=[list(ch1.outfile)],
            gmtfile=ch2.gmtfile,
            configfile=ch2.configfile,
        )
        output = "outdir:dir:{{in.configfile | stem}}.intras-pathwayfeatures"
        lang = config.lang.rscript
        order = 4
        envs = {
            "ncores": config.misc.ncores,
            "fgsea": True,
            "prerank_method": "signal_to_noise",
            "top": 10,
        }
        script = (
            "file://../scripts/scrna_metabolic/MetabolicFeaturesIntraSubsets.R"
        )
        plugin_opts = {
            "report": (
                "file://../reports/scrna_metabolic/"
                "MetabolicFeaturesIntraSubsets.svelte"
            )
        }

    return MetabolicInputs

def main() -> Pipen:
    """Build a pipeline for `pipen run` to run"""
    return Pipen(
        name="scrna-metabolic",
        desc="Metabolic landscape analysis for scRNA-seq data"
    ).set_start(build_processes())
