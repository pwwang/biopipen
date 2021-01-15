"""Network (mathmatics) analysis"""
from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name

pDegree = proc_factory(
    desc='List the degree of nodes, order descendingly.',
    config=Diot(annotate="""
    @name:
        pDegree
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn2}}.degree.txt',
    lang=opts.python,
    args=Diot(
        inopts=Diot(),
        infmt='pair-complete',  # matrix,
        cutoff=0,
    )
)
