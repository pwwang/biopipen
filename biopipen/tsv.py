"""TSV file operations"""

from diot import Diot
from .utils import fs2name
from . import opts, proc_factory

# pylint: disable=invalid-name

pMatrixR = proc_factory(
    desc='Operate a matrix and save the new matrix to file',
    config=Diot(annotate="""
    @name:
        pMatrixR
    @description:
        Operate a matrix and save the new matrix to file.
    @input:
        `infile:file`: The input file containing the matrix
    @output:
        `outfile:file`: The output matrix
    @args:
        `inopts`: The input options for infile:
            - `cnames`: Whether the input file has cnames. Default: True
            - `rnames  `: Whether the input file has rnames. Default: True
            - `delimit`: The delimit. Default: `\t`
            - `skip`: First N lines to skip. Default: `0`
        `params`: Other params for `read.table`. Default: `{"check.names": "FALSE", "quote": ""}`
        `code`: The R code to operating the matrix. (the matrix is read in variable `mat`)
    """),
    lang=opts.Rscript,
    input="infile:file",
    output="outfile:file:{{i.infile | bn}}",
    args=Diot(inopts=Diot(cnames=True, rnames=True, delimit="\t", skip=0),
              params=Diot({
                  "check.names": "FALSE",
                  "quote": ""
              }),
              code=[])
)

pTranspose = proc_factory(
    desc='Transpose a matrix',
    config=Diot(annotate="""
    @name:
        pTranspose
    @description:
        Transpose a matrix
    @input:
        `infile:file`: The input matrix file
    @output:
        `outfile:file`: The Transposed matrix file. Default: `{{i.infile | bn}}`
    @args:
        `inopts`: Input options for input file.
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | bn}}',
    lang=opts.Rscript,
    args=Diot(inopts=Diot(cnames=True, rnames=True))
)

pPaired = proc_factory(
    desc="Subset each input file and make sure they have paired columns.",
    config=Diot(annotate="""
    @name:
        pPaired
    @description:
        Subset each input file and make sure they have paired columns.
    @input:
        `infile1:file`: The first file
        `infile2:file`: The second file
    @outfile:
        `outfile1:file`: The paired file for infile1. Default: `{{i.infile1 | fn}}.paired{{i.infile1 | ext}}`
        `outfile2:file`: The paired file for infile2. Default: `{{i.infile2 | fn}}.paired{{i.infile2 | ext}}`
    @args:
        `inopts1`: reading options for input file1
        `inopts2`: reading options for input file2
    """),
    lang=opts.python,
    input='infile1:file, infile2:file',
    output=[
        'outfile1:file:{{i.infile1 | fn}}.paired{{i.infile1 | ext}}',
        'outfile2:file:{{i.infile2 | fn}}.paired{{i.infile2 | ext}}'
    ],
    args=Diot(
        inopts1=Diot(head=True, headCallback=None),
        inopts2=Diot(head=True, headCallback=None),
    )
)

pCbind = proc_factory(
    desc='Cbind the rest of files to the first file.',
    config=Diot(annotate="""
    @input:
        infiles: The input files
    @output:
        outfile: The output matrix
    @args:
        inopts (Diot): The input options for infile
            - Options for read.table.inopts in utils/__init__.r
        fn2cname (str): An R function used to convert file name to column name.
            - It can have 1 or 2 arguments
            - If only 1 argument is given, it is the filename (without path) of each input file.
            - If 2 arguments are given, 1st is the filename and 2nd is the original column names.
        fill: Do `cbind.fill` instead of `cbind`.
            - Set it to `False` if the row names are in the same order
        na: Replacement for missing values.
    """),
    lang=opts.Rscript,
    input='infiles:files',
    output='outfile:file:{{i.infiles | [0] | stem2 }}_etc.cbound.txt',
    args=Diot(inopts=Diot(cnames=True, rnames=True),
              na='NA',
              fn2cname='function(fn) fn',
              fill=True)
)
pTsvCbind = pCbind.copy()

pRbind = proc_factory(
    desc='Rbind the rest of files to the first file.',
    config=Diot(annotate="""
    @name:
        pRbind
    @description:
        Rbind the rest of files to the first file.
    @input:
        `infiles:files`: The input files
    @output:
        `outfile:file`: The output matrix
    @args:
        `inopts`: The input options for infile:
            - `cnames`: Whether the input file has cnames. Default: True
            - or [True, True, False] corresponding to the file order
            - `rnames  `: Whether the input file has rnames. Default: True
            - `delimit`: The delimit. Default: `\t`
            - `skip`: First N lines to skip. Default: `0`
        `params`: Other params for `read.table`. Default: `{"check.names": "FALSE", "quote": ""}`
        `na`: Replacement for missing values. Default: `NA`
        `fn2rname`: The function (r) used to convert file name to row name.
        `fill`: Do `rbind.fill` instead of `rbind`. Default: `True`
            - Set it to `False` if the row names are in the same order
        `na`: Replacement for missing values. Default: `NA`
    """),
    input='infiles:files',
    output='outfile:file:{{i.infiles[0] | bn}}',
    lang=opts.Rscript,
    args=Diot(
        inopts=Diot(
            # or [True, True, False] corresponding to the file order
            cnames=True,
            rnames=True,
            delimit="\t",
            skip=0
        ),
        params=Diot({"check.names": "FALSE", "quote": ""}),
        na='NA',
        fn2rname='function(fn) fn',
        fill=True,
    ),
    envs=Diot(fs2name=fs2name)
)

pTsvRbind = pRbind.copy()

pCsplit = proc_factory(
    desc='Split the columns of input file into different files.',
    config=Diot(annotate="""
    @input:
        infile: The input file
    @output:
        outdir: The directory containing the output column files
    @args:
        inopts (Diot): The input options for infile:
            - Options for read.table.inopts in utils/__init__.r
        by (str|int): How should the columns be split. Could be:
            - An integer as the chuck size
            - An R function to group the columns. If colnames are not available, indexes will be used.
                - Accept one column name each time. So if-else may be used lots of times
    """),
    input='infile:file',
    output='outdir:dir:{{i.infile | fn}}.csplits',
    lang=opts.Rscript,
    args=Diot(inopts=Diot(cnames=True, rnames=True), by=1)
)

pRsplit = proc_factory(
    desc='Rbind the rest of files to the first file.',
    config=Diot(annotate="""
                @description:
                    Split a matrix by rows and save them into files.
                @input:
                    infile: The input file
                @output:
                    outdir: The directory containing the output row files
                @args:
                    inopts (Diot): The input options for infile
                    size: The chunk size (how many rows to split into one file)
                """),
    input='infile:file',
    output='outdir:dir:{{i.infile | fn}}.rsplits',
    lang=opts.Rscript,
    args=Diot(
        inopts=Diot(cnames=True, rnames=True, delimit="\t", skip=0),
        size=1,
    )
)

pTsv = proc_factory(
    config=Diot(annotate="""
    @name:
        pTsv
    @description:
        Read, Transform, filter a TSV file.
    @input:
        `infile:file`: The input file
    @output:
        `outfile:file`: The output file
    @args:
        `inopts`: The input options for infile:
        `outopts`: The output options for outfile
        `row`: A row function to transform/filter the row. Argument is an instance of `TsvRecord`
        `helper`: A helper function for `args.ops`
    """),
    desc='Read, Transform, filter a TSV file.',
    input="infile:file",
    output="outfile:file:{{i.infile | fn}}.tsv",
    lang=opts.python,
    args=Diot(helper='',
              row=None,
              inopts=Diot(delimit='\t',
                          comment='#',
                          skip=0,
                          cnames=True),
              outopts=Diot(delimit='\t', cnames=True))
)

pTsvDump = proc_factory(
    desc="Dump certain column from a TSV file",
    config=Diot(annotate="""
    @input:
        infile: The input file
    @output:
        outfile: The output file
    @args:
        inopts (Diot): The input options for infile
        col (int|str): The index of the column (0-based) or the column name if inopts.cnames is True
        unique (bool): Whether output the unique values only
    """),
    lang=opts.python,
    input='infile:file',
    output='outfile:file:{{i.infile | stem}}.dumped.txt',
    args=Diot(inopts=Diot(cnames=True),
              col=0,
              unique=True)
)

pTsvColFilter = proc_factory(
    desc='Filter a tsv file by columns',
    config=Diot(annotate="""
    @input:
        infile : The input file
        colfile: The file with columns, one per line, or a list of columns separated by comma.
            - If this is provided, `args.cols` will be ignored.
            - If columns are from a file, then they should be column names instead of indexes.
    @output:
        outfile: The output file
    @args:
        inopts: The options for reading input file. Default: `Diot(cnames = True)`
        keep  : Whether to keep in `args.cols` or to discard
        cols  : The columns used to filter. Could be:
            - column names
            - indices(0-based)
            - a string of lambda function to select column names
            - a file containing the column names, one per line.
    """),
    lang=opts.python,
    input='infile:file, colfile:var',
    output='outfile:file:{{i.infile | bn}}',
    args=Diot(
        inopts=Diot(cnames=True),
        keep=True,
        cols=None,
    )
)
pTsvColSelect = pTsvColFilter.copy()

pTsvAggregate = proc_factory(
    config=Diot(annotate="""
        @input:
            infile: The input file
        @output:
            outfile: The output file
                - With columns `args.on` and aggregated results from `args.aggrs`
                - If `args.on` is a function, then the calculated term will be add to the 1st column.
        @args:
            inopts (Diot): The options to read the input file
            on (str|int): Aggregate according to which column, Default: `0`
                - It also could column name if `args.inopts = True`
                - The input file has to sorted by this column
                - Or a (lambda) function to calculate the term to aggregate on.
            helper (str): Raw codes to give some help for `args.aggrs`
            origin (str): Whether keep or drop the original columns.
                - `drop`: Drop the original columns
                - `keep`: Keep all original records
                - `keep1`: Keep the first original record
                - `keep0`: alias of `keep1`
                - `keep-1`: Keep the last original record
            sorted (bool): Whether the input file has been sorted by `args.on`
            aggrs (Diot): The aggregation methods. Required.
                - It's a `Diot` with the keys for aggregated results
                - If `args.inopts.cnames = True` then the keys will be output as column names, otherwise ignored
                - You can also combine the aggregation results.
                    - For example: `{"sum,mean": lambda rs: [sum(r.value for r in rs), sum(r.value for r in rs)/float(len(rs))]}`
                - We have some built-in aggregation methods:
                    - `args.aggrs.Sum = "$sum:3"`: Get the sum of the 4th column
                    - `args.aggrs.First = "$first:ID"`: Get the first record of the ID column
                    - `args.aggrs.Last = "$last:ID"`: Get the last record of the ID column
                    - `args.aggrs.Mean = "$mean:Height`: Get the mean of column "Height"
                    - `args.aggrs.Median = "$median:1"`: Get the median of the 2nd column
                    - `args.aggrs.Min = "$min:1"`: Get the min of the 2nd column
                    - `args.aggrs.Max = "$max:1"`: Get the max of the 2nd column
                    - `args.aggrs.Max2 = "$max:2"`: Get the max of the 3rd column
                    - `args.aggrs["Max,Max2"] = "$max:1,2"`: Get the max of the 2nd and 3rd column, respectively
                    - `args.aggrs.CombinedP = "$fisher:1"`: Get the combined pvalues for 1st column using fisher'method (`scipy.stats.combine_pvalues`)
        """),
    desc='Aggregate on columns with a set of records',
    input='infile:file',
    output='outfile:file:{{i.infile | fn2}}.aggr.txt',
    lang=opts.python,
    args=Diot(
        inopts=Diot(cnames=True),
        on=0,  # which column
        aggrs=Diot(),
        origin='drop',  # keep
        sorted=True,
        helper=''
    )
)

pTsvHeader = proc_factory(
    desc='Get the header of a tsv file.',
    config=Diot(annotate="""
    @name:
        pTsvHeader
    @description:
        Get the header of a TSV file
    @input:
        `infile:file`: The input file
    @output:
        `outfile:file`: The output file, Default: `{{i.infile | fn2}}.header.txt`
    @args:
        `inopts`: The options to read input file. Default: `Diot(cnames = True)`
        `filter`: The filter for the header. Default: `None`
            - `None`: no filter
            - `lambda cnames: ...` A callback to manipulate colnames.
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn2}}.header.txt',
    lang=opts.python,
    args=Diot(
        inopts=Diot(cnames=True),
        filter=None,
    )
)

pTsvReplaceHeader = proc_factory(
    desc="Replace the header of a tsv file.",
    config=Diot(annotate="""
    @name:
        pTsvReplaceHeader
    @description:
        Replace the header of a TSV file
    @input:
        `infile:file`: The input file
        `hfile:file`:  The file containing the headers, one per line.
    @output:
        `outfile:file`: The output file, Default: `{{i.infile | bn}}`
    @args:
        `inopts`: The options to read input file, Default: `Diot(cnames = True)`
        `cnames`: The column names or callback, Default: `None`
            - `None`: use the header in `i.hfile`
            - `<list/str/file>`: the header to use if `i.hfile` is not provided
            - `lambda cnames, bname: ...`: The callback to modify header in `i.hfile` if provided, otherwise modify the original header with the basename of `i.infile`.
    """),
    lang=opts.python,
    input='infile:file, hfile:file',
    output='outfile:file:{{i.infile | bn}}',
    args=Diot(inopts=dict(cnames=True), cnames=None),
)

pTsvJoin = proc_factory(
    config=Diot(
        echo_jobs=0,
        annotate="""
        @name:
            pTsvJoin
        @description:
            Read files simultaneously.
            NOTE: only one file allows multiple lines with same value to compare, and that file should be the first one. For example:
            ```
            File1:
            1	1
            1	2
            1	3
            File2:
            1	1
            2	2
            3	3
            ```
            If you compare the first column, File1 has to put at the begining for input.
        @input:
            infiles:files: The input files
        @output:
            outfile:file: The output file
        @args:
            inopts: The input options for infile, each can be a list if different for `infiles`:
                - skip   : First N lines to skip. Default: `0`
                - delimit: The delimit. Default: `\t`
                - comment: The comment line mark. Default: `#`
                - cnames : Whether input file has head. Default: `True`
            outopts: The output options:
                - delimit: The delimit. Default: `\t`
                - cnames : Whether to output the head? Default: `False`
                - headCallback: a string callback to organize headers.
                    - E.g: `lambda cnames: ','.join(cnames)`
                    - The argument `cnames` will be the joint cnames from all readers
                    - Could be `True` to use `\t` to assemble the headers
            debug: Save debug information in stderr file. Default: `False`
            match: The match function.
                - Return -1 if matched, otherwise return index the of record that should be read next.
                - E.g: `lambda r1, r2: -1 if r1.id == r2.id else 0 if r1.id < r2.id else 1`
                - By default, will use the first column to compare
            do: The do function. Global vaiable `fout` is available to write results to output file.
                - E.g: `lambda writer, r1, r2: writer.write(r1)`
            helper: Some helper codes.
        """
    ),
    desc='Read files simultaneously.',
    input='infiles:files',
    output='outfile:file:{{i.infiles[0] | fn}}.etc.joined.txt',
    lang=opts.python,
    args=Diot(
        inopts=Diot(delimit='\t', skip=0, comment='#', cnames=True),
        outopts=Diot(delimit='\t', cnames=False),
        debug=False,
        match=None,
        do=None,
        helper='',
    )
)

pTsvSql = proc_factory(
    desc='Query tsv file using SQL',
    config=Diot(annotate="""
    @name:
        pTsvSql
    @description:
        Query tsv file using SQL. (see: http://harelba.github.io/q/examples.html)
    @input:
        `infile:file` : The input tsv file
        `sqlfile:file`: The file containing the SQLs. If provided, `args.sql` will be ignored.
    @output:
        `outfile:file`: The output file
    @args:
        `sql`: If SQL to execute. Use `-` for table name
        `inopts`: Options for input file.
            - `cnames`: Input file has header? Default: `True`
            - `delimit`: The delimit of input file. Default: `\t`
            - `encoding`: Encoding of input file. Default: `UTF-8`
            - `gz`: Whether input file is gzipped. Default: `auto` (detected from file extension)
        `outopts`: Output options.
            - `cnames`: Inherited from `args.inopts`
            - `delimit`: Inherited from `args.inopts`
            - `encoding`: Inherited from `args.inopts`
    @requires:
        q:
          desc: A command line tool that allows direct execution of SQL-like queries on CSVs/TSVs
          url: http://harelba.github.io/q/index.html
          version: 1.7.1
          validate: "{{args.harelba_q}} -v"
          install: |
            wget https://cdn.rawgit.com/harelba/q/1.7.1/bin/q -O /tmp/q;
            chmod +x /tmp/q;
            sed -i 's@#!/usr/bin/env python@#!/usr/bin/env python2@' /tmp/q;
            install /tmp/q "{{args.harelba_q | ?.__contains__: '/' | =:_ | !:'$bindir$/'}}"
    """),
    input='infile:file, sqlfile:file',
    output='outfile:file:{{i.infile | fn}}.bysql{{i.infile | ext}}',
    lang=opts.python,
    args=Diot(
        harelba_q=opts.harelba_q,
        sql='',
        inopts=Diot(cnames=True,
                    delimit="\t",
                    encoding='UTF-8',
                    gz='auto'),
        outopts=Diot(cnames=None, delimit=None, encoding=None)
    )
)

pTsvSample = proc_factory(
    desc='Sample records from a TSV file.',
    lang=opts.python,
    config=Diot(annotate="""
    @name:
        pTsvSample
    @description:
        Sample records from a TSV file
    @input:
        infile: The input file
    @output:
        outfile: The output file
    @args:
        inopts   : input options, only skip available
        n        : how many records to sample`
        arsample : sample program by Alex Reynolds
        replace  : Whether sample with replacement or not
        keeporder: Keep the order of the sampled records as it's in input file
        seed     : The seed
        params   : Other params for arsample
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn2}}.sampled.txt',
    args=Diot(
        inopts=Diot(),
        n=10,
        arsample=opts.arsample,
        replace=False,
        keeporder=False,
        seed=0,
        params=Diot(),
    )
)

pTsvMerge = proc_factory(
    desc='Merge files by rows.',
    config=Diot(annotate="""
    @name:
        pTsvMerge
    @description:
        Merge files in the input directory
    @input:
        `indir:file`: The input directory
    @output:
        `outfile:file`: The output file
    @args:
        `inopts`: The options for input file. Default: `Diot(skip = 0, comment = '#', delimit = '\t')`
        `outopts`: The options for output file. Default: `Diot()`
    """),
    input="infiles:files",
    output="outfile:file:{{i.infiles | fs2name}}",
    lang=opts.python,
    args=Diot(
        inopts=Diot(skip=0, comment='#', delimit='\t'),
        outopts=Diot(),
        # IOError: [Errno 24] Too many open files
        maxopen=100,
    ),
    envs=Diot(fs2name=fs2name)
)

pMergeRows = proc_factory(
    desc='Merge repeated rows.',
    config=Diot(annotate="""
    @name:
        pMergeRows
    @description:
        Merge repeated rows
    @input:
        `infile:file`: The input file (has to be sorted by the repeated columns)
    @output:
        `outfile:file`: The output file. Default: `{{i.infile | bn}}`
    @args:
        `inopts`: The options for input file.
            - defaults: skip: 0, comment: #, delimit '\\t'
        `outopts`: The options for output file. Defaults:
            - head: False (not output head line)
            - headPrefix: `#` (The prefix for head line)
            - headDelimit: `\\t` (The delimit for head line)
            - headTransform: `None` (The callback for head line)
            - delimit: `\\t` (The delimit for output line)
        `match`: The function to return a value to decide whether the row is repeated, argument is a `TsvRecord`.
        `do`   : The merge function in python, argument is a list of `TsvRecord`s or `list`s if `args.inopts.ftype` is `nometa`
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | bn}}',
    lang=opts.python,
    args=Diot(
        inopts=Diot(skip=0, comment='#', delimit='\t'),
        outopts=Diot(head=False,
                     headPrefix='',
                     headDelimit='\t',
                     headTransform=None,
                     delimit='\t'),
        match=None,
        do=None
    )
)

pTsvSplit = proc_factory(
    desc='Split a tsv file by rows.',
    config=Diot(annotate="""
                @description:
                    Split a TSV by rows.
                    If you have `Too many open files` error, using `pRsplit`
                    instead, which will load the whole data into memory.
                @input:
                    infile: The input TSV file
                @output:
                    outdir: The output directory of split files
                @args:
                    inopts (Diot): Options for TsvReader to read input file
                    outopts (Diot): Options for TsvWriter to write output file
                    by (str|int): Split by what?
                        - Size as an integer
                        - Column as "col:column_name" or "col:column_index" (0-based)
                        - A string of lambda function using TsvRecord as argument, returns a tag of split file.
                """),
    lang=opts.python,
    input="infile:file",
    output="outdir:dir:{{i.infile | fn}}.splits",
    args=Diot(inopts=Diot(cnames=True),
              outopts=Diot(header=True),
              by=None)
)

pTsvs2Xlsx = proc_factory(
    desc='Save tsv files to xlsx sheets.',
    config=Diot(annotate="""
    @name:
        pTsvs2Xlsx
    @description:
        Save tsv files to xlsx sheets.
    @input:
        `infiles:files`: The input tsv files
    @output:
        `outfile:file`: The output xlsx file
    @args:
        `fn2sheet`: How to convert filename(without extension) to sheet name
    @requires:
        openpyxl:
          desc    : A Python library to read/write Excel 2010 xlsx/xlsm files
          url     : https://openpyxl.readthedocs.io/en/stable/
          version : 3.0.0
          install : "{{proc.lang}} -m pip install openpyxl"
          validate: "{{proc.lang}} -c 'import openpyxl; from packaging import version; \
              version.parse(openpyxl.__version__) >= version.parse(\\"3.0.0\\")'"
    """),
    input='infiles:files',
    output='outfile:file:{{i.infiles | fs2name}}.xlsx',
    lang=opts.python,
    args=Diot(fn2sheet=None),
    envs=Diot(fs2name=fs2name)
)

pTsv2Xlsx = proc_factory(
    desc='Convert a single TSV file to XLSX file',
    config=Diot(annotate="""
    @input:
        infile: the input TSV file
    @output:
        outfile: the output XLSX file
    @requires:
        openpyxl:
          desc    : A Python library to read/write Excel 2010 xlsx/xlsm files
          url     : https://openpyxl.readthedocs.io/en/stable/
          version : 3.0.0
          install : "{{proc.lang}} -m pip install openpyxl"
          validate: "{{proc.lang}} -c 'import openpyxl; from packaging import version; \
              version.parse(openpyxl.__version__) >= version.parse(\\"3.0.0\\")'"
    """),
    lang=opts.python,
    input='infile:file',
    output='outfile:file:{{i.infile | stem}}.xlsx'
)

pTsvsFromXlsx = proc_factory(
    desc="Convert each XLSX sheet to TSV file.",
    config=Diot(annotate="""
    @input:
        infile: The input XLSX file
    @output:
        outdir: The output directory of converted TSV files
    @requires:
        openpyxl:
          desc    : A Python library to read/write Excel 2010 xlsx/xlsm files
          url     : https://openpyxl.readthedocs.io/en/stable/
          version : 3.0.0
          install : "{{proc.lang}} -m pip install openpyxl"
          validate: "{{proc.lang}} -c 'import openpyxl; from packaging import version; \
              version.parse(openpyxl.__version__) >= version.parse(\\"3.0.0\\")'"
    """),
    lang=opts.python,
    input="infile:file",
    output="outdir:dir:{{i.infile | stem}}.tsvs",
)

pTsvFromXlsx = proc_factory(
    config=Diot(annotate="""
    @input:
        infile: The input XLSX file
    @output:
        outfile: The output converted TSV file
    @requires:
        openpyxl:
        desc    : A Python library to read/write Excel 2010 xlsx/xlsm files
        url     : https://openpyxl.readthedocs.io/en/stable/
        version : 3.0.0
        install : "{{proc.lang}} -m pip install openpyxl"
        validate: "{{proc.lang}} -c 'import openpyxl; from packaging import version; \
            version.parse(openpyxl.__version__) >= version.parse(\\"3.0.0\\")'"
    """),
    desc="Convert active XLSX sheet to TSV file.",
    lang=opts.python,
    input="infile:file",
    output="outfile:file:{{i.infile | stem}}.tsv",
)
