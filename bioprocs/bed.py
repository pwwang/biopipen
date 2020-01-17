"""Processes for BED files"""
from diot import Diot
from . import params, proc_factory
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
        bedtools=params.bedtools.value,
        bedops=params.bedops_sort.value,
        mem='8G',
        by='coord',
        unique=True,
        params=Diot(),
        chrorder=None,
        tmpdir=params.tmpdir.value,
    ),
    lang=params.python.value,
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
    lang=params.python.value,
    args=Diot(
        liftover=params.liftover.value,
        lochain=params.lochain.value,
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
    lang=params.python.value,
    input='infile:file',
    output='outfile:file:{{i.infile | stem}}.bed',
    args=Diot(bedcols=Diot(), keepattrs=True, outhead='#'),
)

pBedFromGff = pGff2Bed.copy()
