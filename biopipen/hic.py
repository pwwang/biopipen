"""
Some long range interaction processes, not necessarily Hi-C
"""
from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name

pPartners = proc_factory(
    desc='Find the interaction partners of the regions in input file.',
    config=Diot(annotate="""
    @name:
        pPartners
    @description:
        Find the interaction partners of the regions in input file.
    @input:
        `regfile:file`: The region file for the regions to find partners for.
        `intfile:file`: The interaction file
    @output:
        `outfile:file`: The regions with partners.
    @args:
        `regtype`: The type of region file. Default: `auto` (tell from file extension)
            - Could also be `bed` or `bedx`
        `inttype`: The type of interaction file. Default: `auto`
            - Could also be `bedpe`, `chiapet.tool`, `hiclib` and `bed12`
    """),
    input="regfile:file, intfile:file",
    output="outfile:file:{{i.regfile | fn2}}-{{i.intfile | fn2}}.partners.bedx",
    lang=opts.python,
    args=Diot(
        regopts=Diot(ftype="auto"),  # bed,   bedx
        intopts=Diot(ftype="auto")  # bedpe, chiapet.tool, hiclib, bed12
    )
)
