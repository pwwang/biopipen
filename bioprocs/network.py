"""Network (mathmatics) analysis"""
from pyppl import Proc
from diot import Diot
from . import params, proc_factory

pDegree = proc_factory(
	desc = 'List the degree of nodes, order descendingly.',
	config = Diot(annotate = """
	@name:
		pDegree
	"""))
pDegree.input       = 'infile:file'
pDegree.output      = 'outfile:file:{{i.infile | fn2}}.degree.txt'
pDegree.args.inopts = Diot()
pDegree.args.infmt  = 'pair-complete' # matrix
pDegree.args.cutoff = 0
pDegree.lang        = params.python.value
