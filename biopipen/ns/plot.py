"""Plotting data"""

import warnings

from ..core.proc import Proc
from ..core.config import config

warnings.warn(
    "The `biopipen.ns.plot` module is deprecated and will be removed in the future. "
    "Please use `biopipen.ns.misc.Plot` process instead.",
    DeprecationWarning,
)


class VennDiagram(Proc):
    """Plot Venn diagram

    Needs `ggVennDiagram`

    Input:
        infile: The input file for data
            If `envs.intype` is raw, it should be a data frame with row names
            as categories and only column as elements separated by comma (`,`)
            If it is `computed`, it should be a data frame with row names
            the elements and columns the categories. The data should be binary
            indicator (`0, 1`) indicating whether the elements are present
            in the categories.

    Output:
        outfile: The output figure file

    Envs:
        inopts: The options for `read.table()` to read `in.infile`
        intype: `raw` or `computed`. See `in.infile`
        devpars: The parameters for `png()`
        args: Additional arguments for `ggVennDiagram()`
        ggs: Additional ggplot expression to adjust the plot
    """

    input = "infile:file"
    output = "outfile:file:{{in.infile | stem}}.venn.png"
    lang = config.lang.rscript
    envs = {
        "inopts": {"row.names": -1, "header": False},
        "intype": "raw",
        "devpars": {"res": 100, "width": 800, "height": 600},
        "args": {},
        "ggs": None,
    }
    script = "file://../scripts/plot/VennDiagram.R"


class Heatmap(Proc):
    """Plot heatmaps using `ComplexHeatmap`

    Examples:
        >>> pipen run plot Heatmap \
        >>> --in.infile data.txt \
        >>> --in.annofiles anno.txt \
        >>> --envs.args.row_names_gp 'r:fontsize5' \
        >>> --envs.args.column_names_gp 'r:fontsize5' \
        >>> --envs.args.clustering_distance_rows pearson \
        >>> --envs.args.clustering_distance_columns pearson \
        >>> --envs.args.show_row_names false \
        >>> --envs.args.row_split 3 \
        >>> --args.devpars.width 5000 \
        >>> --args.devpars.height 5000 \
        >>> --args.draw.merge_legends \
        >>> --envs.args.heatmap_legend_param.title AUC \
        >>> --envs.args.row_dend_reorder \
        >>> --envs.args.column_dend_reorder \
        >>> --envs.args.top_annotation \
        >>>   'r:HeatmapAnnotation( \
        >>>       Mutation = as.matrix(annos[,(length(groups)+1):ncol(annos)]) \
        >>>   )' \
        >>> --envs.args.right_annotation \
        >>>   'r:rowAnnotation( \
        >>>       AUC = anno_boxplot(as.matrix(data), outline = F) \
        >>>   )' \
        >>> --args.globals \
        >>>   'fontsize8 = gpar(fontsize = 12); \
        >>>    fontsize5 = gpar(fontsize = 8); \
        >>>    groups = c  ("Group1", "Group2", "Group3")' \
        >>> --args.seed 8525

    Input:
        infile: The data matrix file
        annofiles: The files for annotation data

    Output:
        outfile: The heatmap plot
        outdir: Other data of the heatmap
            Including RDS file of the heatmap, row clusters and col clusters.

    Envs:
        inopts: Options for `read.table()` to read `in.infile`
        anopts: Options for `read.table()` to read `in.annofiles`
        draw: Options for `ComplexHeatmap::draw()`
        args: Arguments for `ComplexHeatmap::Heatmap()`
        devpars: The parameters for device.
        seed: The seed
        globals: Some globals for the expression in `args` to be evaluated

    Requires:
        bioconductor-complexheatmap:
            - check: {{proc.lang}} <(echo "library(ComplexHeatmap)")
    """
    input = "infile:file, annofiles:files"
    output = [
        'outfile:file:{{in.infile | stem0 | append: ".heatmap"}}/'
        '{{in.infile | stem0 | append: ".heatmap"}}.png',
        'outdir:dir:{{in.infile | stem0 | append: ".heatmap"}}',
    ]
    lang = config.lang.rscript
    envs = {
        "inopts": {"header": True, "row.names": -1},
        "anopts": {"header": True, "row.names": -1},
        "draw": {},
        "devpars": {},
        "args": {"heatmap_legend_param": {}},
        "seed": None,
        "globals": "",
    }
    script = "file://../scripts/plot/Heatmap.R"


class ROC(Proc):
    """Plot ROC curve using [`plotROC`](https://cran.r-project.org/web/packages/plotROC/vignettes/examples.html).

    Input:
        infile: The input file for data, tab-separated.
            The first column should be ids of the records (this is optional if `envs.noids` is True).
            The second column should be the labels of the records (1 for positive, 0 for negative).
            If they are not binary, you can specify the positive label by `envs.pos_label`.
            From the third column, it should be the scores of the different models.

    Output:
        outfile: The output figure file

    Envs:
        noids: Whether the input file has ids (first column) or not.
        pos_label: The positive label.
        ci: Whether to use `geom_rocci()` instead of `geom_roc()`.
        devpars: The parameters for `png()`
        args: Additional arguments for `geom_roc()` or `geom_rocci()` if `envs.ci` is True.
        style_roc: Arguments for `style_roc()`
    """  # noqa: E501
    input = "infile:file"
    output = "outfile:file:{{in.infile | stem}}.roc.png"
    lang = config.lang.rscript
    envs = {
        "noids": False,
        "pos_label": 1,
        "ci": False,
        "devpars": {"res": 100, "width": 750, "height": 600},
        "args": {"labels": False},
        "style_roc": {},
        "show_auc": True,
    }
    script = "file://../scripts/plot/ROC.R"


class Manhattan(Proc):
    """Plot Manhattan plot.

    Using the [`ggmanh`](https://bioconductor.org/packages/devel/bioc/vignettes/ggmanh/inst/doc/ggmanh.html) package.
    Requires `ggmanh` v1.9.6 or later.

    Input:
        infile: The input file for data
            It should contain at least three columns, the chromosome, the position
            and the p-value of the SNPs.
            Header is required.

    Output:
        outfile: The output figure file

    Envs:
        chrom_col: The column for chromosome
            An integer (1-based) or a string indicating the column name.
        pos_col: The column for position
            An integer (1-based) or a string indicating the column name.
        pval_col: The column for p-value
            An integer (1-based) or a string indicating the column name.
        label_col: The column for label.
            Once specified, the significant SNPs will be labeled on the plot.
        devpars (ns): The parameters for `png()`
            - res (type=int): The resolution
            - width (type=int): The width
            - height (type=int): The height
        title: The title of the plot
        ylabel: The y-axis label
        rescale (flag): Whether to rescale the p-values
        rescale_ratio_threshold (type=float): Threshold of that triggers the rescale
        signif (auto): A single value or a list of values to indicate the significance levels
            Multiple values should be also separated by comma (`,`).
            The minimum value will be used as the cutoff to determine if the SNPs are significant.
        hicolors (auto): The colors for significant and non-significant SNPs
            If a single color is given, the non-significant SNPs will be in grey.
            Set it to None to disable the highlighting.
        thin_n (type=int): Number of max points per horizontal partitions of the plot.
            `0` or `None` to disable thinning.
        thin_bins (type=int): Number of bins to partition the data.
        zoom (auto): Chromosomes to zoom in
            Each chromosome should be separated by comma (`,`) or in a list. Single chromosome is also accepted.
            Ranges are also accepted, see `envs.chroms`.
            Each chromosome will be saved in a separate file.
        zoom_devpars (ns): The parameters for the zoomed plot
            - width (type=int): The width
            - height (type=int): The height, inherited from `devpars` by default
            - res (type=int): The resolution, inherited from `devpars` by default
        chroms (auto): The chromosomes and order to plot
            A hyphen (`-`) can be used to indicate a range.
            For example `chr1-22,chrX,chrY,chrM` will plot all autosomes, X, Y and M.
            if `auto`, only the chromosomes in the data will be plotted in the order
            they appear in the data.
        args (ns): Additional arguments for `manhattan_plot()`.
            See <https://rdrr.io/github/leejs-abv/ggmanh/man/manhattan_plot.html>.
            Note that `-` will be replaced by `.` in the argument names.
            - <more>: Additional arguments for `manhattan_plot()`
    """  # noqa: E501
    input = "infile:file"
    output = "outfile:file:{{in.infile | stem0}}.manhattan.png"
    lang = config.lang.rscript
    envs = {
        "chrom_col": 1,
        "pos_col": 2,
        "pval_col": 3,
        "label_col": None,
        "devpars": {"res": 100, "width": 1000, "height": 500},
        "zoom_devpars": {"width": 500, "height": None, "res": None},
        "title": None,
        "ylabel": "-log10(p-value)",
        "rescale": True,
        "rescale_ratio_threshold": 5,
        "signif": [5e-8, 1e-5],
        "hicolors": None,
        "thin_n": None,
        "thin_bins": 200,
        "zoom": None,
        "chroms": "auto",
        "args": {},
    }
    script = "file://../scripts/plot/Manhattan.R"


class QQPlot(Proc):
    """Generate QQ-plot or PP-plot using qqplotr.

    See <https://cran.r-project.org/web/packages/qqplotr/vignettes/introduction.html>.

    Input:
        infile: The input file for data
            It should contain at least one column of p-values or the values to be
            plotted. Header is required.
        theorfile: The file for theoretical values (optional)
            This file should contain at least one column of theoretical values.
            The values will be passed to `envs.theor_qfunc` to calculate the theoretical
            quantiles.
            Header is required.

    Output:
        outfile: The output figure file

    Envs:
        val_col: The column for values to be plotted
            An integer (1-based) or a string indicating the column name.
        devpars (ns): The parameters for `png()`
            - res (type=int): The resolution
            - width (type=int): The width
            - height (type=int): The height
        xlabel: The x-axis label
        ylabel: The y-axis label
        title: The title of the plot
        trans: The transformation of the values
            You can use `-log10` to transform the values to `-log10(values)`.
            Otherwise you can a direct R function or a custom R function.
            For example `function(x) -log10(x)`.
        kind (choice): The kind of the plot, `qq` or `pp`
            - qq: QQ-plot
            - pp: PP-plot
        theor_col: The column for theoretical values in `in.theorfile` if provided,
            otherwise in `in.infile`.
            An integer (1-based) or a string indicating the column name.
            If `distribution` of `band`, `line`, or `point` is `custom`, this column
            must be provided.
        theor_trans: The transformation of the theoretical values.
            The `theor_funs` have default functions to take the theoretical values.
            This transformation will be applied to the theoretical values before
            passing to the `theor_funs`.
        theor_funs (ns): The R functions to generate density, quantile and deviates
            of the theoretical distribution base on the theoretical values
            if `distribution` of `band`, `line`, or `point` is `custom`.
            - dcustom: The density function, used by band
            - qcustom: The quantile function, used by point
            - rcustom: The deviates function, used by line
        args (ns): The common arguments for `envs.band`, `envs.line` and `envs.point`.
            - distribution: The distribution of the theoretical quantiles
                When `custom` is used, the `envs.theor_col` should be provided and
                `values` will be added to `dparams` automatically.
            - dparams (type=json): The parameters for the distribution
            - <more>: Other shared arguments between `stat_*_band`, `stat_*_line`
                and `stat_*_point`.
        band (ns): The arguments for `stat_qq_band()` or `stat_pp_band()`.
            See <https://rdrr.io/cran/qqplotr/man/stat_qq_band.html> and
            <https://rdrr.io/cran/qqplotr/man/stat_pp_band.html>.
            Set to `None` or `band.disabled` to True to disable the band.
            - disabled (flag): Disable the band
            - distribution: The distribution of the theoretical quantiles
                When `custom` is used, the `envs.theor_col` should be provided and
                `values` will be added to `dparams` automatically.
            - dparams (type=json): The parameters for the distribution
            - <more>: Additional arguments for `stat_qq_band()` or `stat_pp_band()`
        line (ns): The arguments for `stat_qq_line()` or `stat_pp_line()`.
            See <https://rdrr.io/cran/qqplot/man/stat_qq_line.html> and
            <https://rdrr.io/cran/qqplot/man/stat_pp_line.html>.
            Set to `None` or `line.disabled` to True to disable the line.
            - disabled (flag): Disable the line
            - distribution: The distribution of the theoretical quantiles
                When `custom` is used, the `envs.theor_col` should be provided and
                `values` will be added to `dparams` automatically.
            - dparams (type=json): The parameters for the distribution
            - <more>: Additional arguments for `stat_qq_line()` or `stat_pp_line()`
        point (ns): The arguments for `geom_qq_point()` or `geom_pp_point()`.
            See <https://rdrr.io/cran/qqplot/man/stat_qq_point.html> and
            <https://rdrr.io/cran/qqplot/man/stat_pp_point.html>.
            Set to `None` or `point.disabled` to True to disable the point.
            - disabled (flag): Disable the point
            - distribution: The distribution of the theoretical quantiles
                When `custom` is used, the `envs.theor_col` should be provided and
                `values` will be added to `dparams` automatically.
            - dparams (type=json): The parameters for the distribution
            - <more>: Additional arguments for `geom_qq_point()` or `geom_pp_point()`
        ggs (list): Additional ggplot expression to adjust the plot.
    """
    input = "infile:file, theorfile:file"
    output = "outfile:file:{{in.infile | stem}}.{{envs.kind}}.png"
    lang = config.lang.rscript
    envs = {
        "val_col": 1,
        "theor_col": None,
        "theor_trans": None,
        "theor_funs": {
            "dcustom": """
              function(x, values, ...) {
                density(values, from = min(values), to = max(values), n = length(x))$y
              }
            """,
            "qcustom": "function(p, values, ...) {quantile(values, probs = p)}",
            "rcustom": "function(n, values, ...) { sample(values, n, replace = TRUE) }",
        },
        "args": {"distribution": "norm", "dparams": {}},
        "devpars": {"res": 100, "width": 1000, "height": 1000},
        "xlabel": "Theoretical Quantiles",
        "ylabel": "Observed Quantiles",
        "title": "QQ-plot",
        "trans": None,
        "kind": "qq",
        "band": {"disabled": False, "distribution": None, "dparams": None},
        "line": {"disabled": False, "distribution": None, "dparams": None},
        "point": {"disabled": False, "distribution": None, "dparams": None},
        "ggs": None,
    }
    script = "file://../scripts/plot/QQPlot.R"


class Scatter(Proc):
    """Generate scatter plot using ggplot2.

    [`ggpmisc`](https://cran.r-project.org/web/packages/ggpmisc/index.html) is used
    for the stats and labels.
    See also https://cran.r-project.org/web/packages/ggpmisc/vignettes/model-based-annotations.html

    Input:
        infile: The input file for data
            It should contain at least two columns for x and y values.
            Header is required.

    Output:
        outfile: The output figure file

    Envs:
        x_col: The column for x values
            An integer (1-based) or a string indicating the column name.
        y_col: The column for y values
            An integer (1-based) or a string indicating the column name.
        devpars (ns): The parameters for `png()`
            - res (type=int): The resolution
            - width (type=int): The width
            - height (type=int): The height
        args (ns): Additional arguments for `geom_point()`
            See <https://ggplot2.tidyverse.org/reference/geom_point.html>.
            - <more>: Additional arguments for `geom_point()`
        mapping: Extra mapping for all geoms, including `stats`.
            Should be `aes(color = group)` but all these are valid: `color = group` or
            `(color = group)`.
        ggs (list): Additional ggplot expression to adjust the plot.
        formula: The formula for the model
        stats (type=json): The stats to add to the plot.
            A dict with keys available stats in `ggpmisc` (without `stat_`).
            See <https://cran.r-project.org/web/packages/ggpmisc/vignettes/model-based-annotations.html#statistics>.
            The values should be the arguments for the stats.
            If you want a stat to be added multiple times, add a suffix `#x` to the key.
            For example, `poly_line#1` and `poly_line#2` will add two polynomial lines.
    """  # noqa: E501
    input = "infile:file"
    output = "outfile:file:{{in.infile | stem}}.scatter.png"
    lang = config.lang.rscript
    envs = {
        "x_col": 1,
        "y_col": 2,
        "devpars": {"res": 100, "width": 1000, "height": 800},
        "args": {},
        "mapping": None,
        "ggs": [],
        "formula": "y ~ x",
        "stats": {},
    }
    script = "file://../scripts/plot/Scatter.R"
