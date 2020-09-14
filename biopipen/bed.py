"""Processes for BED files"""
from diot import Diot
from . import opts, proc_factory
from .bedtools import * # pylint: disable=wildcard-import,unused-wildcard-import

# pylint: disable=invalid-name

pBedSort = proc_factory(
    desc='Sort bed files.',
    config=Diot(annotate="""
    @name:
        pBedSort
    @description:
        Sort bed files
    @input:
        `infile:file`: The input file
    @output:
        `outfile:file`: The output file
    @args:
        `tool`:     The tool used to sort the file. Default: sort (bedtools, bedops)
        `bedtools`: The path to bedtools. Default: bedtools
        `bedops`:   The path to bedops' sort-bed. Default: sort-bed
        `mem`:      The memory to use. Default: 8G
        `by`:       Sort by coordinates("coord", default) or name("name")
            - Only available when use tool `sort`
        `tmpdir`:   The tmpdir to use. Default: `$TMPDIR`
        `unique`:   Remove the dupliated records? Default: True
        `params`:   Other params for `tool`. Default: {}
        `chrorder`: The chromosome order used to sort. Default: `None`
            - `None`: Sort by natural order (chr1 followed by chr10, instead of chr2)
            - Only available when using `sort` (`args.tool = 'sort'`)
    @requires:
        [`bedtools`](http://bedtools.readthedocs.io/en/latest/index.html)
        [`bedops`](https://github.com/bedops/bedops)
    """),
    input="infile:file",
    output="outfile:file:{{i.infile | bn}}",
    args=Diot(
        tool='sort',
        bedtools=opts.bedtools,
        bedops=opts.bedops_sort,
        mem='8G',
        by='coord',
        unique=True,
        params=Diot(),
        chrorder=None,
        tmpdir=opts.tmpdir,
    ),
    lang=opts.python,
)

pBedSwitchBase = proc_factory(
    desc='Switch BED files between 0-based and 1-based coordinates',
    config=Diot(annotate="""
                @input:
                    infile: The input BED file
                @output:
                    outfile: The output BED file
                @args:
                    inbase (int): The input BED coordinate base
                        - If `args.outbase` is specified as 0, this is implied 1, vice versa.
                    outbase (int): The output BED coordinate base
                        - This is also implied if `args.inbase` is specified
                """),
    input='infile:file',
    output='outfile:file:{{i.infile | bn}}',
    lang=opts.python,
    args=Diot(inbase=None, outbase=None)
)

pBedLiftover = proc_factory(
    desc='Lift over bed files.',
    config=Diot(annotate="""
    @name:
        pBedLiftover
    @description:
        Lift over bed files.
    @input:
        `infile:file`: The input bed file
    @output:
        `outfile:file`: The output file
        `umfile:file` : The unmapped file
    @args:
        `liftover`: The liftover program
        `lochain` : the liftover chain file
    @require:
        `liftover` from UCSC
    """),
    input='infile:file',
    output=['outfile:file:{{i.infile | bn}}',
            'umfile:file:{{i.infile | fn}}.unmapped{{i.infile | ext}}'],
    lang=opts.python,
    args=Diot(
        liftover=opts.liftover,
        lochain=opts.lochain,
        params=Diot(),
    )
)

pGff2Bed = proc_factory(
    desc='Convert GTF/GFF file to BED file',
    config=Diot(annotate="""
    @input:
        infile: The input gtf/gff file
    @output:
        outfile: The converted bed file
    @args:
        bedcols:  Strings of python functions used to convert GTF/GFF records to BED fields.
            - You can define the NAME column here, and extra columns after the 6th column.
            - For example: `args.bedcols = {"NAME": "lambda attrs: rec.CHR + ':' + rec.START"}`
                - `attrs` are the attributes of GFF records, plus CHR, START, END, SCORE and STRAND.
                - See: https://github.com/pwwang/pygff
                - By default, NAME will use `id` in attributes, and then `name`. Otherwise `CHR:START-END` will  be used.
            - You can also add extra columns starting from 7th column of BED file, for example:
                - `args.bedcols = {"CN": "lambda attrs: attrs['CopyNumber']"}`
        keepattrs: Keep the original attributes at last column of the output BED file.
        outhead: Put head to output file or not.
            - Could be prefix to the head.
    """),
    lang=opts.python,
    input='infile:file',
    output='outfile:file:{{i.infile | stem}}.bed',
    args=Diot(bedcols=Diot(), keepattrs=True, outhead='#'),
)

pBedFromGff = pGff2Bed.copy()

pBed2Wiggle = proc_factory(
    desc="Convert BED file to fixedStep Wiggle(bigwig) file",
    config=Diot(annotate="""
                @input:
                    infile: The input BED file
                        - Last column is the signals (separated by comma)
                @output:
                    outfile: The output Wiggle(bigwig) file
                @args:
                    inbase (int): The BED coordinate base of input file
                    step (int): The step of each bin
                    span (int): The span of each bin
                        - Note that the length of of region specified in the input BED file won't be checked if it is exactly Ndata x Msteps
                    addchr (bool): Try to add `chr` prefix to chromosome?
                    bigwig (bool): Output bigwig file? Requires `args.wigtobigwig` and `args.gsize`
                    wigtobigwig (str): Path to ucsc wigToBigWig
                    gsize (str): Path to geome size
                """),
    lang=opts.python,
    input="infile:file",
    output="""outfile:file:{{i.infile | stem}}.{{args.bigwig | ?
                                                             | =:'bw'
                                                             | !:'wig'}}""",
    args=Diot(inbase=0, span=1, step=1, addchr=True,
              gsize=opts.gsize,
              bigwig=True, wigtobigwig=opts.wigtobigwig)
)

pBedExpand2 = proc_factory(
    desc="Expand a bed file with signals at last columns",
    config=Diot(annotate="""
                @description:
                    Expand a bed file with signals at last columns.
                    Extending the `bedtools expand`, which does not make windows
                    of the regions but just splits the signals

                    Bedtools expand behavior:
                    ```
                    $ cat test.txt
                    chr1	10	20	10,20
                    chr1	40	50	40,50

                    $ bedtools expand test.txt -c 5
                    chr1	10	20	10
                    chr1	10	20	20
                    chr1	40	50	40
                    chr1	40	50	50
                    ```

                    This tool's behavior with window size 5:
                    ```
                    chr1	10	15	1,2	10
                    chr1	15	20	1,2	20
                    chr1	40	45	4,5	40
                    chr1	45	50	4,5	50
                    ```
                @input:
                    infile: The input bed file with signals at last column
                        - Signals must be separated by comma
                @output:
                    outfile: The output bed file with both regions and signals
                        being split
                @args:
                    base (int): The coordinate base of input bed file
                    binsize (int): The size of the bins
                        - The length of each region must be divided by `binsize`
                    datacols (str|list): The indexes of the columns that have
                                         the data
                        - 1-based
                        - Number of data must be divided by the number of bins
                          in that region
                    name (str): The name of the output bins
                        - `None`: don't output name
                        - `src`: Use the source name if provided otherwise
                          `src<src_index>` will be used.
                        - `srcbin`: Use the source name plus the bin index
                          (`src_bin<bin_index>`). if src is not provided, using
                          `src<src_index>`
                    origcols (str): What about other columns in the input file
                        - `keep`: Keep and duplicate them for each bin
                        - `drop`: Drop them
                    aggrs (str|list): How to aggregate the data if there are more
                                      than 1 for the bin
                        - `raw`: comma-separated raw values
                        - `mean`: Take the mean (sum/binsize)
                        - `sum`: Take the sum
                        - `min`: Take the min
                        - `max`: Take the max
                        - `count`: Count the number of data in that bin
                        - Arithmetic function will turn the data into a float
                          number
                        - A string of lambda function to aggregate. Argument is
                          the list of data in that bin.
                        - If there are more than 1 datacols, you can specify
                          multiple aggrs in a list for each datacol.
                """),
    lang=opts.python,
    input='infile:file',
    output='outfile:file:{{i.infile | stem}}.expanded.bed',
    args=Diot(base=0, binsize=None, datacols=None,
              name='srcbin', origcols='keep', aggrs='raw')
)
