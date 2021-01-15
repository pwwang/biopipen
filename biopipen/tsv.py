from .core.proc import Proc
from .core.args import args
from .core.defaults import SCRIPT_DIR

MODULE = 'tsv'

class TsvReplaceHeader(Proc):
    """Replace the header of a TSV file"""
    name = 'TsvReplaceHeader'

class TsvColSelect(Proc):
    """Filter/Select columns from a TSV file

    Input:
        infile: The input file
        colfile: The columns to select or a file with columns, one per line.

    Output:
        outfile: The output file

    Args:
        inopts: The input options to read the infile
        keep: Whether keep the provided columns or discard (keep other columns)
        cols: The columns to select, will be overriden by `in.colfile`
    """
    lang = args.python
    input_keys = 'infile:file, colfile'
    output = 'outfile:file:{{in.infile | bn}}'
    script = f'file://{SCRIPT_DIR}/{MODULE}/TsvColSelect.py'
    args = {
        'inopts': {'cnames': True},
        'keep': True,
        'cols': None
    }

class TsvDataFrameR(Proc):
    """Operate the TSV file as a dataframe in R

    dplyr and tidyr operations are supported.

    Input:
        infile: The input file containing the matrix

    Output:
        outfile: The output matrix

    Args:
        inopts: The input options for infile:
            - `cnames`: Whether the input file has cnames. Default: True
            - `rnames  `: Whether the input file has rnames. Default: True
            - `delimit`: The delimit. Default: `\t`
            - `skip`: First N lines to skip. Default: `0`
        code: The R code to operating the dataframe.
            the dataframe is read in variable `df`
    """
    lang = args.Rscript
    input_keys = "infile:file"
    output = "outfile:file:{{in.infile | bn}}"
    script = f'file://{SCRIPT_DIR}/{MODULE}/TsvDataFrameR.R'
    args = dict(
        inopts=dict(cnames=True, rnames=True, delimit="\t", skip=0),
        code=[]
    )

class TsvCBind(Proc):
    """Cbind the rest of files to the first file

    Input:
        infiles: The input files

    Output:
        outfile: The output matrix

    Args:
        inopts: The input options for infile
            - Options for read.table.inopts in utils/__init__.R
        fn2cname: An R function used to convert file name to column name.
            - It can have 1 or 2 arguments
            - If only 1 argument is given, it is the filename (without path) of each input file.
            - If 2 arguments are given, 1st is the filename and 2nd is the original column names.
        fill: Do `cbind.fill` instead of `cbind`.
            - Set it to `False` if the row names are in the same order
        na: Replacement for missing values.
    """
    lang = args.Rscript
    input_keys = 'infiles:files'
    output = 'outfile:file:{{in.infiles | getitem: 0 | stem2 }}_etc.cbound.txt'
    script = f'file://{SCRIPT_DIR}/{MODULE}/TsvCBind.R'
    args = dict(
        inopts=dict(cnames=True, rnames=True),
        na='NA',
        fn2cname='function(fn) fn',
        fill=True
    )

class TsvRBind(Proc):
    """Rbind the rest of files to the first file

    Input:
        infiles: The input files

    Output:
        outfile: The output matrix

    Args:
        inopts: The input options for infile:
            - `cnames`: Whether the input file has cnames.
            - `rnames  `: Whether the input file has rnames.
            - `delimit`: The delimit.
            - `skip`: First N lines to skip.
        na: Replacement for missing values.
        fn2rname: The function (r) used to convert file name to row name.
        fill: Do `rbind.fill` instead of `rbind`. Default: `True`
            - Set it to `False` if the row names are in the same order
    """
    lang = args.Rscript
    input_keys = 'infiles:files'
    output = 'outfile:file:{{in.infiles | getitem: 0 | stem2 }}_etc.rbound.txt'
    script = f'file://{SCRIPT_DIR}/{MODULE}/TsvRBind.R'
    args = dict(
        inopts=dict(cnames=True, rnames=True),
        na='NA',
        fn2rname='function(fn) fn',
        fill=True
    )