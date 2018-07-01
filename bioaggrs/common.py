from pyppl import Aggr
from bioprocs.common import pSort
from bioprocs.tsv import pSimRead


"""
@name:
	aSimRead2
@description:
	Read 2 files simultaneously
@depends:
	pSort1[*] \
                pSimRead[!]
	pSort2[*] /
"""
aSimRead2 = Aggr(
	pSort.copy(id = 'pSort1'),
	pSort.copy(id = 'pSort2'),
	pSimRead,
	depends = False
)
aSimRead2.starts = aSimRead2.pSort1, aSimRead2.pSort2
aSimRead2.ends   = aSimRead2.pSimRead
# depends 
aSimRead2.pSimRead.depends = aSimRead2.pSort1, aSimRead2.pSort2
# delegates
aSimRead2.delegate('args.params1', 'pSort1', 'args.params')
aSimRead2.delegate('args.params2', 'pSort2', 'args.params')
aSimRead2.delegate('args.inopts1', 'pSort1', 'args.inopts')
aSimRead2.delegate('args.inopts2', 'pSort2', 'args.inopts')
aSimRead2.delegate('args', 'pSimRead', 'args')
aSimRead2.delegate('args1', 'pSort1', 'args')
aSimRead2.delegate('args2', 'pSort2', 'args')
# input
aSimRead2.pSimRead.input = lambda ch1, ch2: [list(row) for row in ch1.cbind(ch2)]