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
aSimRead2.delegate('args.params1', 'pSort1', 'args.params')
aSimRead2.delegate('args.params2', 'pSort2', 'args.params')
aSimRead2.delegate('args.delimit1', 'pSort1', 'args.delimit')
aSimRead2.delegate('args.delimit2', 'pSort2', 'args.delimit')
aSimRead2.delegate('args.skip1', 'pSort1', 'args.skip')
aSimRead2.delegate('args.skip2', 'pSort2', 'args.skip')
aSimRead2.delegate('args.*', 'pSimRead')
# input
aSimRead2.pSimRead.input = lambda ch1, ch2: [
	list(r) for r in ch1.repRow(max(ch1.length(), ch2.length())).cbind(ch2)
]