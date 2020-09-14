"""Tabix utilities"""
from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name

pTabix = proc_factory(
    desc='Use tabix to extract information.',
    config=Diot(annotate="""
    @name:
        pTabix
    @description:
        Use tabix to extract information.
    @input:
        `infile`: a local or remote file
        `region`: a region or a file containing regions
    @output:
        `outfile:file`: The information extracted from the input file
    @args:
        `tabix`: The path to `tabix`
        `params`: Other params for `tabix`
    """),
    input="infile, region",
    output="outfile:file:{{i.infile | fn | fn}}-{{job.index}}{{\
            i.infile | ?.endswith: '.gz' | =:[:-3] | $ext}}",
    lang=opts.python,
    args=Diot(
        tabix=opts.tabix,
        params=Diot(h=True),
    )
)

pTabixIndex = proc_factory(
    desc='Generate tabix index file',
    config=Diot(annotate="""
    @name:
        pTabixIndex
    @description:
        Generate tabix index file.
    @input:
        `infile:file`: the input file
            - Could be bgzipped.
    @output:
        `outfile:file`: The bgzipped file
        `outidx:file`:  The tabix index file
    @args:
        `tabix`: The path to `tabix`
        `params`: Other params for `tabix`
    """),
    lang=opts.python,
    input="infile:file",
    output=[
        "outfile:file:{{i.infile | bn}}{% if args.gz %}.gz{% endif %}",
        "outidx:file:{{i.infile  | bn}}{% if args.gz %}.gz{% endif %}.tbi"
    ],
    args=Diot(
        gz=True,
        tabix=opts.tabix,
        params=Diot(),
    )
)
