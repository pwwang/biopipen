from pyppl import Aggr
from bioprocs.common import pSort
from bioprocs.matrix import pSimRead


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
	pSort.copy(newid = 'pSort1'),
	pSort.copy(newid = 'pSort2'),
	pSimRead,
	depends = False
)
aSimRead2.starts = aSimRead2.pSort1, aSimRead2.pSort2
aSimRead2.ends   = aSimRead2.pSimRead
# defaults
aSimRead2.pSort1.args.noeline = True
aSimRead2.pSort2.args.noeline = True
# depends 
aSimRead2.pSimRead.depends = aSimRead2.pSort1, aSimRead2.pSort2
# delegates
aSimRead2.delegate('args.skip', 'pSimRead')
aSimRead2.delegate('args.delimit', 'pSimRead')
aSimRead2.delegate('args.gzip', 'pSimRead')
aSimRead2.delegate('args.match', 'pSimRead')
aSimRead2.delegate('args.do', 'pSimRead')
aSimRead2.delegate('args.usehead', 'pSimRead')
# input
aSimRead2.pSimRead.input = lambda ch1, ch2: [list(r) for r in ch1.cbind(ch2)]