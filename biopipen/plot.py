"""Generate plots using given data"""

from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name

pPlot = proc_factory(
    desc='Generate plot using ggplot2',
    config=Diot(annotate="""
    @name:
        pPlot
    @description:
        Use ggplot2 to generate plots
    @input:
        `infile:file`: The input data file
    @output:
        `outfile:file`: The output file
    @args:
        `cnames` : Whether the input file has colnames. Default: True
        `rnames` : Whether the input file has rownames. Default: False
        `aes`    : The default aes. Default: {'x':1, 'y':2} (corresponding to colnames)
        `helper` : Some helper codes to generate `params` and `ggs`
        `devpars`: The device parameters. Default: `Diot(res = 300, height = 2000, width = 2000)`
        `ggs`    : The extra ggplot elements.
    @requires:
        ggplot2:
          desc: A system for declaratively creating graphics, based on The Grammar of Graphics.
          url: https://ggplot2.tidyverse.org/index.html
          install: "{{proc.lang}} <(echo install.packages('ggplot2'))"
          validate: "{{proc.lang}} <(echo require('ggplot2'))"
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}.png',
    lang=opts.Rscript,
    args=Diot(
        cnames=True,
        rnames=False,
        aes=Diot(x=1, y=2),  # only allow x, y,
        helper='',
        devpars=Diot(res=300, height=2000, width=2000),
        params=Diot(),
        ggs=Diot(),
    )
)

pScatter = proc_factory(
    desc='Generate scatter plot.',
    config=Diot(annotate="""
    @name:
        pScatter
    @description:
        Use ggplot2 geom_point to generate plots
    @infile:
        `infile:file`: The input data file
    @outfile:
        `outfile:file`: The output file
    @args:
        `cnames` : Whether the input file has colnames. Default: True
        `rnames` : Whether the input file has rownames. Default: False
        `x`      : The x aes. Default: 1 (corresponding to colnames)
        `y`      : The y aes. Default: 2 (corresponding to colnames)
        `helper` : Some helper codes to generate `params` and `ggs`
        `devpars`: The device parameters. Default: `Diot(res = 300, height = 2000, width = 2000)`
        `params` : The extra params for `geom_point`
        `ggs`    : The extra ggplot elements.
    @requires:
        ggplot2:
          desc: A system for declaratively creating graphics, based on The Grammar of Graphics.
          url: https://ggplot2.tidyverse.org/index.html
          install: "{{proc.lang}} <(echo install.packages('ggplot2'))"
          validate: "{{proc.lang}} <(echo require('ggplot2'))"
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}.scatter.png',
    lang=opts.Rscript,
    args=Diot(
        cnames=True,
        rnames=False,
        x=1,
        y=2,
        helper='',
        devpars=Diot(res=300, height=2000, width=2000),
        params=Diot(),
        ggs=Diot(),
    )
)

pHisto = proc_factory(
    desc='Generate histogram.',
    config=Diot(annotate="""
    @name:
        pHisto
    @description:
        Use ggplot2 geom_histogram to generate histograms
    @input:
        `infile:file`: The input data file
    @output:
        `outfile:file`: The output file
    @args:
        `inopts` : Options to read the input file. Default: `Diot(cnames = True, rnames = False)`
        `x`      : The x aes. Default: 1 (corresponding to colnames)
        `helper` : Some helper codes to generate `params` and `ggs`
        `devpars`: The device parameters. Default: `Diot(res = 300, height = 2000, width = 2000)`
        `params` : The extra params for `geom_histogram`
        `ggs`    : The extra ggplot elements.
    @requires:
        ggplot2:
          desc: A system for declaratively creating graphics, based on The Grammar of Graphics.
          url: https://ggplot2.tidyverse.org/index.html
          install: "{{proc.lang}} <(echo install.packages('ggplot2'))"
          validate: "{{proc.lang}} <(echo require('ggplot2'))"
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}.histo.png',
    lang=opts.Rscript,
    args=Diot(
        inopts=Diot(cnames=True, rnames=False),
        x=1,
        helper='',
        devpars=Diot(res=300, height=2000, width=2000),
        params=Diot(),
        ggs=Diot(),
    )
)

pDensity = proc_factory(
    desc='Generate density plot.',
    config=Diot(annotate="""
    @name:
        pDensity
    @description:
        Use ggplot2 geom_density to generate density plot
    @input:
        `infile:file`: The input data file
    @output:
        `outfile:file`: The output file
    @args:
        `inopts` : Options to read the input file. Default: `Diot(cnames = True, rnames = False)`
        `x`      : The x aes. Default: 1 (corresponding to colnames)
        `helper` : Some helper codes to generate `params` and `ggs`
        `devpars`: The device parameters. Default: `Diot(res = 300, height = 2000, width = 2000)`
        `params` : The extra params for `geom_density`
        `ggs`    : The extra ggplot elements.
    @requires:
        ggplot2:
          desc: A system for declaratively creating graphics, based on The Grammar of Graphics.
          url: https://ggplot2.tidyverse.org/index.html
          install: "{{proc.lang}} <(echo install.packages('ggplot2'))"
          validate: "{{proc.lang}} <(echo require('ggplot2'))"
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}.density.png',
    lang=opts.Rscript,
    args=Diot(
        inopts=Diot(cnames=True, rnames=False),
        x=1,
        helper='',
        devpars=Diot(res=300, height=2000, width=2000),
        params=Diot(),
        ggs=Diot(),
    )
)

pFreqpoly = proc_factory(
    desc='Generate frequency polygon plot.',
    config=Diot(annotate="""
    @name:
        pFreqpoly
    @description:
        Use ggplot2 geom_freqpoly to generate frequency polygon plot.
    @infile:
        `infile:file`: The input data file
    @outfile:
        `outfile:file`: The output file
    @args:
        `cnames` : Whether the input file has colnames. Default: True
        `rnames` : Whether the input file has rownames. Default: False
        `x`      : The x aes. Default: 1 (corresponding to colnames)
        `helper` : Some helper codes to generate `params` and `ggs`
        `devpars`: The device parameters. Default: `Diot(res = 300, height = 2000, width = 2000)`
        `params` : The extra params for `geom_point`
        `ggs`    : The extra ggplot elements.
    @requires:
        ggplot2:
          desc: A system for declaratively creating graphics, based on The Grammar of Graphics.
          url: https://ggplot2.tidyverse.org/index.html
          install: "{{proc.lang}} <(echo install.packages('ggplot2'))"
          validate: "{{proc.lang}} <(echo require('ggplot2'))"
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}.freqpoly.png',
    lang=opts.Rscript,
    args=Diot(
        cnames=True,
        rnames=False,
        x=1,
        helper='',
        devpars=Diot(res=300, height=2000, width=2000),
        params=Diot(),
        ggs=Diot(),
    )
)

pBoxplot = proc_factory(
    desc='Generate boxplot plot.',
    config=Diot(annotate="""
    @name:
        pBoxplot
    @description:
        Generate box plot
    @input:
        `infile:file`: The data file
    @output:
        `outpng:file`: The output figure
    @args:
        `inopts` :   Input options to read the input file
            - `cnames` :   Whether the input file has header. Default: `True`
            - `rnames` :   Whether the input file has row names. Default: `False`
            - `delimit`:   The seperator. Defualt: `\\t`
        `x`      :   The `ind` (index) column. Only for `args.stacked = True`. Default: `2`
        `y`      :   The `values` column. Only for `args.stacked = True`. Default: `1`
        `helper` :   Some raw codes to help to construct the matrix and arguments.
        `stacked`:   Whether the input file is stacked
            - Stacked file looks like:
            ```
            values	ind
            1.1	col1
            1.2	col1
            ...
            .8	col2
            .9	col2
            ...
            3.2	col3
            ...
            ```
            - Unstacked file looks like:
            ```
            col1	col2	col3
            1.1	.8	3.2
            1.2	.9	2.2
            ```
        `params`:    Other parameters for `geom_boxplot`, default: `Diot()`
        `ggs`   :    Extra ggplot2 statements
    @requires:
        ggplot2:
          desc: A system for declaratively creating graphics, based on The Grammar of Graphics.
          url: https://ggplot2.tidyverse.org/index.html
          install: "{{proc.lang}} <(echo install.packages('ggplot2'))"
          validate: "{{proc.lang}} <(echo require('ggplot2'))"
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}.boxplot.png',
    lang=opts.Rscript,
    args=Diot(
        inopts=Diot(cnames=True, rnames=False, delimit='\t'),
        x=2,
        y=1,
        helper='',
        stacked=False,
        devpars=Diot(res=300, height=2000, width=2000),
        params=Diot(),
        ggs=Diot(),
    )
)

pBar = proc_factory(
    desc='Generate bar/col plot.',
    config=Diot(annotate="""
    @name:
        pBar
    @description:
        Generate bar/col plot
    @input:
        `infile:file`: The data file
    @output:
        `outpng:file`: The output figure
    @args:
        `inopts` :   Input options to read the input file
            - `cnames` :   Whether the input file has header. Default: `True`
            - `rnames` :   Whether the input file has row names. Default: `False`
            - `delimit`:   The seperator. Defualt: `\\t`
        `x`      :   The `ind` (index) column. Only for `args.stacked = True`. Default: `2`
        `y`      :   The `values` column. Only for `args.stacked = True`. Default: `1`
        `helper` :   Some raw codes to help to construct the matrix and arguments.
        `stacked`:   Whether the input file is stacked
            - see `pBoxplot.args.stacked`
        `params`:    Other parameters for `geom_bar`, default: `Diot()`
        `ggs`   :    Extra ggplot2 statements
    @requires:
        ggplot2:
          desc: A system for declaratively creating graphics, based on The Grammar of Graphics.
          url: https://ggplot2.tidyverse.org/index.html
          install: "{{proc.lang}} <(echo install.packages('ggplot2'))"
          validate: "{{proc.lang}} <(echo require('ggplot2'))"
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}.bar.png',
    lang=opts.Rscript,
    args=Diot(
        inopts=Diot(cnames=True, rnames=False, delimit='\t'),
        x=2,
        y=1,
        helper='',
        stacked=False,
        devpars=Diot(res=300, height=2000, width=2000),
        params=Diot(),
        ggs=Diot(),
    )
)

pHeatmap = proc_factory(
    desc='Plot heatmaps.',
    config=Diot(annotate="""
    @name:
        pHeatmap
    @description:
        Plot heatmaps.
    @input:
        `infile:file`: The input matrix file
    @output:
        `outfile:file`: The heatmap
    @args:
        `ggs`: The ggplot items for heatmap
        `devpars`: The parameters for device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`
        `dendro`: The parameters for control of the dendrogram. Default: `{'dendro': True}`
            - `dendro`: `True`: plot dendros for both rows and cols; `col`: only plot dendro for cols; `row`: only plot dendro for rows
            - `rows`: The rownames to subset the rows and control the order of rows. Must a list. Only works when not plotting dendro for rows.
            - `cols`: The colnames to subset the cols and control the order of cols. Must a list. Only works when not plotting dendro for cols.
        `header`: The input file has header? Default: True
        `rownames`: The input file has rownames? Default: 1
        `rows`: Row selector
            - `all`: All rows
            - `top:N`: Top N rows (original data ordered in descending order). N defaults to 100
            - `bottom:N`: Bottom N rows. N defaults to 100
            - `both:N`: Top N rows and bottom N rows. N defaults to 50
            - `random:N`: Random N rows. N defaults to 50
            - `random-both:N`: Random N rows from top part and N rows from bottom part. N defaults to 50
        `cols`: Col selector (see `rows`).
    @requires:
        ggplot2:
          desc: A system for declaratively creating graphics, based on The Grammar of Graphics.
          url: https://ggplot2.tidyverse.org/index.html
          install: "{{proc.lang}} <(echo install.packages('ggplot2'))"
          validate: "{{proc.lang}} <(echo require('ggplot2'))"
    """),
    input="infile:file",
    output="outfile:file:{{i.infile | fn}}.heatmap.png",
    lang=opts.Rscript,
    args=Diot(
        ggs=Diot(),
        devpars=Diot(res=300, height=2000, width=2000),
        params=Diot(dendro=True),
        inopts=Diot(rnames=True, cnames=True),
        helper='',
    )
)

pHeatmap2 = proc_factory(
    desc='Plot heatmaps using R package ComplexHeatmap.',
    config=Diot(annotate="""
    @name:
        pHeatmap2
    @description:
        Plot heatmaps using R package ComplexHeatmap. Example:
        ```bash
        bioprocs plot.pHeatmap2
        -i.infile MMPT.txt
        -i.annofiles:l:o PatientAnno.txt
        -args.params.row_names_gp 'r:fontsize5'
        -args.params.column_names_gp 'r:fontsize5'
        -args.params.clustering_distance_rows pearson
        -args.params.clustering_distance_columns pearson
        -args.params.show_row_names false
        -args.params.row_split 3
        -args.devpars.width 5000
        -args.devpars.height 5000
        -args.draw.merge_legends
        -args.params.heatmap_legend_param.title AUC
        -args.params.row_dend_reorder
        -args.params.column_dend_reorder
        -args.params.top_annotation ' \\
            r:HeatmapAnnotation(Mutation = as.matrix(
                                    annos[,(length(groups)+1):ncol(annos)]
                                ),
                                Group = as.matrix(annos[,groups]),
                                col = list(Mutation = c(`0`="grey",
                                                        `1`="lightgreen",
                                                        `2`="green",
                                                        `3`="darkgreen")),
                                annotation_name_gp = fontsize8,
                                show_legend = c(Group=F))'
        -args.params.right_annotation 'r:rowAnnotation(AUC = anno_boxplot(as.matrix(data),
                                                       outline = F))'
        -args.helper 'fontsize8 = gpar(fontsize = 12);
                      fontsize5 = gpar(fontsize = 8);
                      groups = c  ("Group1", "Group2", "Group3")'
        -args.seed 8525
        ```
    @input:
        `infile:file`: The input data file for the main heatmap.
        `annofiles:files`: The annotation files.
            - For now, they should share the same `args.anopts`
    @output:
        `outfile:file`: The plot.
        `outdir:dir`: The output directory including the output plot and other information.
    @args:
        `devpars` : The parameters for device.
        `draw`    : The parameters for `ComplexHeatmap::draw`
        `params`  : Other parameters for `ComplexHeatmap::Heatmap`
        `anopts`  : The options to read the annotation files.
        `inopts`  : The options to read the input files.
        `seed`    : The seed. Default: `None`
        `helper`  : The raw R codes to help defining some R variables or functions. Default: ''
        `saveinfo`: Save other information include the heatmap object and the clusters. Default: `True`
    @requires:
        ComplexHeatmap:
          desc: efficient to visualize associations between different sources of data sets and reveal potential patterns.
          url: https://github.com/jokergoo/ComplexHeatmap
          install: "{{proc.lang}} <(echo devtools::install_github('jokergoo/ComplexHeatmap', ref = 'c269eb4'))"
          validate: "{{proc.lang}} <(echo require('ComplexHeatmap'))"
    """),
    lang=opts.Rscript,
    input="infile:file, annofiles:files",
    output=[
        "outfile:file:{{i.infile | fn2}}.heatmap/{{i.infile|fn2}}.heatmap.png",
        "outdir:dir:{{i.infile | fn2}}.heatmap"
    ],
    args=Diot(
        devpars=Diot(res=300, height=2000, width=2000),
        draw=Diot(),
        params=Diot(heatmap_legend_param=Diot()),
        anopts=Diot(rnames=True, cnames=True),
        inopts=Diot(rnames=True, cnames=True),
        seed=None,
        saveinfo=True,
        helper='',
    )
)

pScatterCompare = proc_factory(
    desc='Plot scatter compare plots.',
    config=Diot(annotate="""
    @input:
        infile: The input file
            - If `args.stacked` is `True`, it will be like (`args.rnames` should be `False`):
              ```
              Item  Value   Group
              rs1   1       A
              rs2   2       A
              ...
              rs1   7       B
              rs2   9       B
              ```
            - Otherwise:
              ```
              [ROWNAME  ]A  B
              [rs1      ]1  7
              [rs2      ]2  9
              ...
              ```
    @output:
        outfile: The output plot
    @args:
        ggs (Diot): Extra expressions for ggplot. Note if geom_point is included, original geom_point will be ignored.
        params (Diot): Other parameters for `geom_point`
        devpars (Diot): The parameters for plot device.
        x (int|str): The (unstacked) column used to be x-axis.
        tsform (str): An R function to transform the data for each group
        corr (str): Method to calculate correlation to show. `None` to hide.
            - Could be one of `pearson`, `spearman` or `kendall`
            - If `args.tsform` is provided, correlation will be calculated based on transformed values.
        line (str): The line to separate the performance.
            - Note that the coefficient and `x` should be connected with `*`
            - For example: `y = 2*x + 1`
            - Points with `y > 2*x + 1` and `y < 2*x + 1` will showed in different colors.
        stacked (bool): Whether the input data is stacked.
        inopts (Diot): Options to read input file.
    """),
    input="infile:file",
    output="outfile:file:{{i.infile | fn}}.scattercomp.png",
    lang=opts.Rscript,
    args=Diot(
        ggs=Diot(),
        params=Diot(),
        devpars=Diot(res=300, height=2000, width=2000),
        x=1,
        tsform=None,
        corr="pearson",
        line="y=x",
        stacked=False,
        inopts=Diot(cnames=True, rnames=True)
    )
)

pROC = proc_factory(
    desc='ROC curves and AUCs.',
    config=Diot(annotate="""
    @input:
        infile: The input matrix file.
            - Col0: rownames if `args.inopts.rnames` is True
            - Col1: label (0, 1 class)
            - Col2: prediction values from model1
            - [Col3: prediction values from model2]
            - [...]
    @output:
        outdir: The output directory
    @args:
        `inopts`: The options for input file. Default: `Diot(rnames = True, cnames = True)`
        `params`: The parameters for `plot.roc` from `utils/plot.r`
        `ggs`   : Additaional ggplot terms. Default:
          ```python
          Diot({
              'style_roc': {},
              # show legend at bottom right corner
              'theme#auc': {'legend.position': [1, 0], 'legend.justification': [1, 0]}
          })
          ```
        `devpars`: The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`
    """),
    lang=opts.Rscript,
    input='infile:file',
    output=[
        'outfile:file:{{i.infile | stem}}.roc/{{i.infile | stem}}.roc.png',
        'outdir:dir:{{i.infile | stem}}.roc'
    ],
    args=Diot(
        inopts=Diot(rnames=True, cnames=True),
        params=Diot(bestCut=True, showAUC=True),
        ggs=Diot({
            # show legend at bottom right corner
            'theme#auc': {
                'legend.position': [1, 0],
                'legend.justification': [1, 0],
                'panel.border':
                'r:element_rect(colour="black", fill=NA, size=.2)'
            }
        }),
        devpars=Diot(res=300, height=2000, width=2000)
    )
)

pVenn = proc_factory(
    desc='Venn plots.',
    config=Diot(annotate="""
    @input:
        infile: The input matrix, could be two formats:
            - `args.intype == "raw"`:
                ```
                category1	category2	category3
                e1	e2	e2
                e2	e3	e4
                ... ...
                ```
                - `args.intype == "computed"`:
                ```
                    category1	category2	category3
                [e1]	0	1	1
                [e2]	0	0	1
                ...
                [eN]	1	0	0
                ```
        metafile: The metadata file for each category for upset plot.
            - format:
                ```
                    col1	col2	...	colN
                category1	x1	y1	...	z1
                category2	x2	y2	...	z2
                ...	...
                categoryN	xN	yN	...	zN
                ```
    @output:
        outfile: The plot
    @args:
        tool (str)    : Which tools to use (`venn`, `upsetr` or `auto`).
            - auto : venn if `n <= 3` otherwise upsetr
        intype (str)  : Type of input file. See `i.infile`.
        inopts (Diot) : options to read the input file.
            - rnames: Whether input file has rownames.
        params (Diot) : Other params for `venn.diagram` or `upset`.
        devpars (Diot): The parameters for plot device.
    """),
    lang=opts.Rscript,
    input="infile:file, metafile:file",
    output="outfile:file:{{i.infile | fn}}.venn.png",
    args=Diot(
        tool='auto',  # upsetr or auto: < = 3 venn, else upsetr,
        inopts=Diot(rnames=False, cnames=True),
        intype='raw',  # computed,
        params=Diot(),
        devpars=Diot(res=300, height=2000, width=2000),
    )
)

pVenn2 = proc_factory(
    desc='Venn plots using individual input files',
    config=Diot(annotate="""
    @name:
        pVenn2
    @description:
        Venn plots using individual input files, each of which contains the elements of the category.
    @input:
        `infiles:files`: The input files, each one is a category containing the elements.
            - If it has column name, then it will be used as category name, otherwise
            - the filename (without extension) will be used.
        `metafile:file`: The metadata file for each category for upset plot.
            - format:
            ```
                col1	col2	...	colN
            category1	x1	y1	...	z1
            category2	x2	y2	...	z2
            ...	...
            categoryN	xN	yN	...	zN
            ```
    @output:
        `outfile:file`: The output file, Default: `{{i.infiles | fs2name}}.venn.png`
    @args:
        `tool`    : Which tools to use. Default: auto (venn, upsetr, auto(n<=3: venn, otherwise upsetr))
        `inopts`  : options to read the input file. Default: `Diot(rnames = False, cnames = False)`
        `params`  : Other params for `venn.diagram` or `upset`. Default: `Diot()`
        `devpars` : The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`
        cats (str): An R function to tranform the categories using the filename and cnames if `args.inopts.cnames == True`
    """),
    input="infiles:files, metafile:file",
    output="outfile:file:{{i.infiles[0] | stem}}_etc.venn.png",
    lang=opts.Rscript,
    args=Diot(
        tool='auto',
        inopts=Diot(rnames=False, cnames=False),
        params=Diot(),
        cats=None,
        devpars=Diot(res=300, height=2000, width=2000),
    )
)

pPie = proc_factory(
    desc='Pie chart.',
    config=Diot(annotate="""
    @name:
        pPie
    @description:
        Plot piechart
    @input:
        `infile:file`: The input file. Could be either:
            - Direct numbers of each category (`args.inopts.rnames = True`).
            ```
            Group1	50
            Group2	50
            ```
            - Transposed number of each category (`args.inopts.cnames = True`)
            ```
            Group1	Group2
            50	50
            ```
            - Presence of each items in the category (`args.inopts.cnames = True`)
            ```
                Group1	Group2
            Item1	1	0
            Item2	0	1
            ...
            ```
            - Or a stacked version:
            ```
            Item1	Group1
            Item2	Group2
            ```
    @output:
        `outfile:file`: the output plot
    @args:
        intype (str): The type of the input file. `direct`, `presence` or `stacked`.
            - See `i.infile`
        inopts (Diot): The options for reading the input file
        ggs (Diot): Extra expressions for ggplot.
        devpars (Diot): The parameters for plot device.
        params (Diot): The parameters for `ggpie` for `ggpubr`
    """),
    input="infile:file",
    output="outfile:file:{{i.infile | fn}}.pie.png",
    lang=opts.Rscript,
    args=Diot(
        inopts=Diot(delimit="\t", rnames=False),
        devpars=Diot(res=300, height=2000, width=2000),
        intype="direct",
        params=Diot(),
        ggs=Diot(),
    )
)

pManhattan2 = proc_factory(
    desc="Manhattan plot.",
    config=Diot(annotate="""
                @input:
                    infile: The input file with columns:
                        - Snp name
                        - Chromsome
                        - Position
                        - Pvalue
                    hifile: The file with highlight data:
                        - Snp name
                        - Labels (optional)
                @output:
                    outfile: The output plot file
                @args:
                    hifile (str): The file with highlight data
                        - ignored if `i.hifile` provided
                    devpars (Diot): The parameters for plot device.
                    params (Diot): Parameters for `CMplot`
                        - See `?CMplot` in R for more options.
                """),
    input="infile:file, hifile:file",
    output="outfile:file:{{i.infile | stem2}}.manhattan.jpg",
    lang=opts.Rscript,
    args=Diot(
        inopts=Diot(cnames=True, rnames=False),
        devpars=Diot(res=300, height=2000, width=2000),
        hifile=None,
        params=Diot()
    )
)

pManhattan = proc_factory(
    desc='Manhattan plot.',
    config=Diot(annotate="""
    @name:
        pManhattan
    @description:
        Manhattan plot.
    @input:
        `infile:file`: The input file. First 6 columns should be BED6, and then column:
            - 7: The raw pvalue.
            - 8: The x-axis labels for the records. Optional. If not provided, chromosomes will be used.
            - This file has to be sorted by coordinates.
            - For example:
                ```
                chr19	45604163	45604163	rs7255060	0	+	3.238E-03	+200K
                chr19	45595277	45595277	rs10417101	0	+	3.870E-03	+200K
                chr19	45394336	45394336	rs71352238	0	+	6.440E-03	-50K
                chr19	45615857	45615857	rs6509194	0	+	1.298E-02	+250K
                chr19	45594170	45594170	rs3178166	0	+	3.617E-02	+200K
                chr19	45361574	45361574	rs3852856	0	+	2.070E-02	-100K
                chr19	45220205	45220205	rs57090948	0	+	4.384E-02	-200K
                chr19	45396219	45396219	rs157582	0	+	9.472E-03	-50K
                chr19	45210634	45210634	rs10421830	0	+	1.375E-02	-250K
                chr19	45228502	45228502	rs10422350	0	+	4.121E-02	-200K
                ```
        `hifile:file`: The file with the record names (one per line) to highlight in the plot.
    @output:
        `outfile:file`: The plot. Default: `{{i.infile | fn}}.manht.png`
    @args:
        `inopts` : Options to read the input file. Default: `Diot(cnames = False, rnames = False)`
        `hilabel`: Show the labels of the highlight points. Default: `True`
        `ggs`    : Extra expressions for ggplot.
        `devpars`: The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`
        `gsize`  : The genome sizes file.
            - If `None`, will infer from the snp coordinates
    """),
    input='infile:file, hifile:file',
    output='outfile:file:{{i.infile | fn}}.manht.png',
    lang=opts.Rscript,
    args=Diot(
        inopts=Diot(cnames=False, rnames=False),
        hilabel=True,
        devpars=Diot(res=300, height=2000, width=2000),
        ggs=Diot(),
        gsize=None  # opts.gsize,
    )
)

pQQ = proc_factory(
    desc='Q-Q plot',
    config=Diot(annotate="""
                @input:
                    infile: The input file with one or two columns.
                        - If only one column is given, x-axis will be theoretical values
                @outfile:
                    outfile: The plot file
                @args:
                    inopts (Diot): The options to read the input file
                    devpars (Diot): The parameters for plot device.
                    ggs (Diot): Extra expression for ggplot object
                    tsform (str): An R function to transform the input data
                    params (Diot): Paramters for `geom_scatter`
                """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}.qq.png',
    lang=opts.Rscript,
    args=Diot(
        inopts=Diot(cnames=True, rnames=False),
        devpars=Diot(res=300, height=2000, width=2000),
        ggs=Diot(),
        tsform=None,
        params=Diot(),
    )
)

pPairs = proc_factory(
    desc='Plot pairs with ggpairs from GGally',
    config=Diot(annotate="""
    @name:
        pPairs
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}.pairs.png',
    lang=opts.Rscript,
    args=Diot(
        inopts=Diot(cnames=True, rnames=True),
        devpars=Diot(res=300, height=400, width=400),  # for each cell,
        ggs=Diot(),  # geom_* not available but theme available,
        # other parameters for ggpairs
        params=Diot(upper=Diot(continuous="density"))
    )
)

pVolcano = proc_factory(
    desc='Do volcano plot.',
    config=Diot(annotate="""
    @input:
        infile: The input file to plot the volcano plot
    @output:
        outfile: The output plot
    @args:
        fccut (float): The (log) fold change cutoff
        pcut (float): The pvalue cutoff
        usepval (bool): Whether use pvalue or qvalue.
        hilights (list): The gene to highlight
        devpars (Diot): The device parameters for the plot
        ggs (Diot): The extra ggs for the plot
        fdrpos (str): Where should we put the FDR=xxx label
            - `left` or `right`.
            - Sometimes it may be blocked by the highlights
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}.volcano.png',
    lang=opts.Rscript,
    args=Diot(
        fccut=2,
        pcut=.05,
        usepval=False,
        fdrpos="right",
        hilights=[],
        devpars=Diot(res=300, height=2000, width=2000),
        ggs=Diot(),
    )
)
