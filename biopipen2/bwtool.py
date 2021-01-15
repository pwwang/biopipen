"""BigWig file operations using bwtool"""

from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name

pBWExtract = proc_factory(
    desc='Extract features from bigwig file.',
    config=Diot(annotate="""
        @description:
            Extract features from bigwig file using `bwtool extract`.
            See: https://github.com/CRG-Barcelona/bwtool/wiki/extract
        @input:
            infile: The input bigwig file.
            bedfile: The bed/gff file of regions to extract features from.
        @output:
            outfile: The output features from those regions
        @args:
            bwtool (str): Path to bwtool.
        """),
    lang=opts.python,
    input='infile:file, bedfile:file',
    output='outfile:file:{{i.infile | stem}}.bwtool-extracted.txt',
    args=Diot(bwtool=opts.bwtool)
)

pBWSummary = proc_factory(
    desc='provide some summary stats for each region in a bed file',
    config=Diot(annotate="""
    @description:
        Provide some summary stats for each region in a bed file using `bwtool summary`
        See: https://github.com/CRG-Barcelona/bwtool/wiki/summary
        Only bed file for regions is supported in this process.
    @input:
        infile: The input bigwig file.
        bedfile: The bed/gff file of regions to extract features from.
    @output:
        outfile: The output file.
    @args:
        bwtool (str): The path to bwtool
        params (Diot): Other parameters for `bwtool summary`
    """),
    lang=opts.python,
    input='infile:file, bedfile:file',
    output='outfile:file:{{i.infile | stem}}.bwtool-summarized.txt',
    args=Diot(bwtool=opts.bwtool,
              params=Diot(header=True, keep_bed=True))
)
