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
# delegates
aSimRead2.delegate('args', 'pSimRead')
# depends 
aSimRead2.pSimRead.depends = aSimRead2.pSort1, aSimRead2.pSort2
# input
aSimRead2.pSimRead.input = lambda ch1, ch2: [list(row) for row in ch1.cbind(ch2)]

"""
@name:
	aTsvJoin2
@description:
	An alias of aSimRead2
"""
aTsvJoin2 = aSimRead2.copy()