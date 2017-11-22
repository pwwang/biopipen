import inspect
from pyppl import PyPPL, Proc
from os import path, walk

filedir = None
config  = {'log': {'level': 'basic', 'lvldiff': "-DEBUG"}}

def _getFiledir():
	global filedir
	fn = inspect.getframeinfo(inspect.getouterframes(inspect.currentframe())[2][0])[0]
	fn = path.basename(fn)
	fn = fn[4].lower() + fn[5:-3]
	filedir = path.join(path.realpath(path.dirname(path.dirname(__file__))), 'testfiles', fn)

_getFiledir()

def getfile (name, input = True):
	return path.join(filedir, 'input' if input else 'expect', name)

def procOK(proc, name, test):
	predfile = proc.channel.get()
	exptfile = getfile(name, False)
	if path.isfile(predfile):
		with open(predfile) as fpred, open(exptfile) as fexpt:
			test.assertEqual(fpred.read().splitlines(), fexpt.read().splitlines())
	else:
		for root, _, files in walk(predfile):
			for name in files:
				pfile = path.join(root, name)
				efile = path.join(exptfile, name)
				with open(pfile) as fpred, open(efile) as fexpt:
					test.assertEqual(fpred.read().splitlines(), fexpt.read().splitlines())

def utilTest(input, script, name, tplenvs, test, args = None):
	ends = {
		'.r' : 'Rscript',
		'.py': 'python'
	}
	pTest         = Proc(desc = 'Test utils.', tag=name.split('.')[0])
	pTest.input   = input
	pTest.output  = "outfile:file:outfile"
	pTest.lang    = [ends[k] for k in ends if script.endswith(k)][0]
	pTest.tplenvs = tplenvs
	pTest.args    = {} if not args else args
	pTest.script  = 'file:' + script
	PyPPL(config).start(pTest).run()
	procOK(pTest, name, test)

