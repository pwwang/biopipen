"""Plotting data"""

from ..core.proc import Proc
from ..core.config import config

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
        "devpars": {"res": 100, "width": 1000, "height": 1000},
        "args": {},
        "ggs": None,
    }
    script = "file://../scripts/plot/VennDiagram.R"
