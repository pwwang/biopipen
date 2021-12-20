"""Metabolic landscape analysis for scRNA-seq data

An abstract from https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape

If you have clustering done somewhere else, you could use `replace_clustering()`

Reference:
    [1] Xiao, Zhengtao, Ziwei Dai, and Jason W. Locasale.
    "Metabolic landscape of the tumor microenvironment at
    single cell resolution." Nature communications 10.1 (2019): 1-12.

Start Process:
    MetabolicInputs
"""
from pathlib import Path
from typing import Any, Mapping

from datar.all import tibble, if_else
from pipen import Pipen
from ..core.config import config
from ..core.proc import Proc
from .scrna import SeuratPreparing, SeuratFilter, SeuratClustering, SCImpute


OPTIONS = {
    "clustered": False
}

def build_processes(options: Mapping[str, Any] = None):
    """Build processes for metabolic landscape analysis pipeline"""
    options = options or OPTIONS

    class MetabolicInputs(Proc):
        """Input for the metabolic pathway analysis pipeline for
        scRNA-seq data

        Input:
            metafile: A metafile indicating the metadata
                Two columns are required: `Sample` and `RNADir`
                `Sample` should be the first column with unique identifiers
                for the samples and `RNADir` indicates where the expression
                matrice are.
                Currently only 10X data is supported
            subsetfile: A file with information to subset the cells.
                Each subset of cells will be treated separately.
                See `subsetting` of `in.config`
            groupfile: The group file to group the cells.
                Rows are groups, and columns are cell IDs
                (without sample prefices)
                from samples or a single column `ALL` but with all cell IDs with
                sample prefices. Or it can be served as a config file for
                subsetting the cells.
                See `grouping` of `in.config`
            gmtfile: The GMT file with the metabolic pathways
            config: String of configurations in TOML for the pipeline
                `grouping_name`: The name of the groupings. Default: `Group`
                `grouping`: How the cells should be grouped.
                If `"Input"`, use `in.groupfile` directly to group the cells
                else if `Idents`, use `Idents(srtobj)` as groups. Otherwise
                (`Config`), otherwise it is TOML config with
                `name => condition` pairs. The `condition` will be passed to
                `subset(srtobj, ...)` to get the subset of cells
                `subsetting`: Similar as `grouping`, indicating how to use
                `in.subsetfile` to subset the cells.

        Output:
            metafile: Soft link to `in.metafile`
            subsetfile: Soft link to `in.subsetfile`
            groupfile: Soft link to `in.groupfile`
            gmtfile: Soft link to `in.gmtfile`
            configfile: The config file with `in.config` saved
        """

        input = """
            metafile:file,
            subsetfile:file,
            groupfile:file,
            gmtfile:file,
            config:var
        """
        output = """
            {%- for inkey, inval in in.items() -%}
            {%- if inkey != "config" -%}
                {{- inkey}}:file:{{inval | str | basename -}},
            {%- endif -%}
            {%- endfor -%}
            configfile:file:config.toml
        """
        script = """
            {% addfilter writefile %}
            def writefile(s, outfile):
                with open(outfile, "w") as f:
                    f.write(s)
            {% endaddfilter %}

            {%- for inkey, inval in in.items() -%}
            {%- if inkey != "config" and inval %}
                ln -s {{inval | quote}} {{out[inkey] | quote -}}
            {% elif inkey != "config" and not inval %}
                touch {{out[inkey] | quote -}}
            {%- endif -%}
            {%- endfor %}

            cat > {{out.configfile | quote}} <<EOF
            {{in.config | replace: "`", "\\`"}}
            EOF

        """


    class MetabolicSeuratPreparing(SeuratPreparing):
        if not options["clustered"]:
            requires = MetabolicInputs


    class MetabolicSeuratClustering(SeuratClustering):
        requires = MetabolicSeuratPreparing


    class MetabolicCellSubsets(SeuratFilter):
        if options["clustered"]:
            requires = MetabolicInputs
            input_data = lambda ch: tibble(
                srtobj=ch.metafile,
                filterfile=if_else(
                    [
                        ssfile is None or Path(ssfile).name == "None"
                        for ssfile in ch.subsetfile
                    ],
                    ch.configfile,
                    ch.subsetfile,
                ),
            )
        else:
            requires = MetabolicSeuratClustering, MetabolicInputs
            input_data = lambda ch1, ch2: tibble(
                srtobj=ch1.rdsfile,
                filterfile=if_else(
                    [
                        ssfile is None or Path(ssfile).name == "None"
                        for ssfile in ch2.subsetfile
                    ],
                    ch2.configfile,
                    ch2.subsetfile,
                ),
            )


    class MetabolicCellGroups(Proc):
        """Group cells for metabolic landscape analysis

        Each group of cells will do the enrichment
        against the metabolic pathways
        """

        requires = MetabolicCellSubsets, MetabolicInputs
        input = "srtdir:file, groupfile:file, configfile:file"
        output = "outdir:dir:{{in.srtdir | stem}}"
        lang = config.lang.rscript
        script = "file://../scripts/scrna_metabolic/MetabolicCellGroups.R"


    class MetabolicExprImputation(SCImpute):
        """Impute the dropout values in scRNA-seq data."""

        requires = MetabolicCellGroups
        input = "srtdir:file"
        output = "outdir:dir:imputed"
        lang = config.lang.rscript
        envs = {"ncores": config.misc.ncores, "drop_thre": 0.5}
        script = "file://../scripts/scrna_metabolic/MetabolicExprImputation.R"


    class MetabolicPrepareSCE(Proc):
        """Prepare SingleCellExperiment objects"""

        requires = MetabolicExprImputation, MetabolicCellGroups, MetabolicInputs
        input = "impdir:dir, srtdir:dir, gmtfile:file"
        output = "outfile:file:metabolic.sce.RDS"
        lang = config.lang.rscript
        envs = {"refexon": config.ref.refexon}
        script = "file://../scripts/scrna_metabolic/MetabolicPrepareSCE.R"


    class MetabolicExprNormalization(Proc):
        """Normalize the expression data using deconvolution"""

        requires = MetabolicPrepareSCE, MetabolicInputs
        input = "sceobj:file, configfile:file"
        output = "outfile:file:{{in.sceobj | stem}}.sce.RDS"
        envs = {"dropout": 0.75, "refexon": config.ref.refexon}
        lang = config.lang.rscript
        script = (
            "file://../scripts/scrna_metabolic/MetabolicExprNormalization.R"
        )


    class MetabolicPathwayActivity(Proc):
        """Pathway activities for each group"""

        requires = MetabolicExprNormalization, MetabolicInputs
        input = "sceobj:file, gmtfile:file, configfile:file"
        output = "outdir:dir:{{in.sceobj | stem}}.pathwayactivity"
        envs = {
            "ntimes": 5000,
            "ncores": config.misc.ncores,
            "heatmap_devpars": {"res": 100, "width": 1200, "height": 2000},
            "violin_devpars": {"res": 100, "width": 1000, "height": 550},
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
        """Pathway heterogeneity"""

        requires = MetabolicExprNormalization, MetabolicInputs
        input = "sceobj:file, gmtfile:file, configfile:file"
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
        """Inter-subset metabolic features - Enrichment analysis in details"""

        requires = MetabolicExprNormalization, MetabolicInputs
        input = "sceobj:file, gmtfile:file, configfile:file"
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
        """Intra-subset metabolic features - Enrichment analysis in details"""

        requires = MetabolicExprNormalization, MetabolicInputs
        input = "sceobj:file, gmtfile:file, configfile:file"
        output = "outdir:dir:{{in.sceobj | stem}}.intras-pathwayfeatures"
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

MetabolicLandscape = (
    Pipen(desc="Metabolic landscape analysis for scRNA-seq data")
    .set_start(build_processes())
)
