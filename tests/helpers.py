import inspect, gzip
from pyppl import PyPPL, Proc
from os import path, listdir

filedir = None
config  = {'log': {'level': 'basic', 'lvldiff': "-DEBUG"}}

def _getFiledir():
	global filedir
	fn = inspect.getframeinfo(inspect.getouterframes(inspect.currentframe())[2][0])[0]
	fn = path.basename(fn)
	fn = fn[4].lower() + fn[5:-3]
	filedir = path.join(path.realpath(path.dirname(path.dirname(__file__))), 'testfiles', fn)

_getFiledir()

def getfile (name = '', input = True):
	return path.join(filedir, 'input' if input else 'expect', name)

def getbin (name = ''):
	return path.join(path.dirname(path.dirname(path.dirname(filedir))), 'bin', name)

def procOK(proc, name, test, order = True):
	predfile = proc.channel.get()
	exptfile = getfile(name, False)
	if path.isfile(predfile):
		openfunc = gzip.open if exptfile.endswith('.gz') else open
		with openfunc(predfile) as fpred, openfunc(exptfile) as fexpt:
			if order:
				test.assertEqual(fpred.read().splitlines(), fexpt.read().splitlines())
			else:
				test.assertEqual(sorted(fpred.read().splitlines()), sorted(fexpt.read().splitlines()))
	else:
		for item in listdir(predfile):
			pfile = path.join(predfile, item)
			efile = path.join(exptfile, item)
			if not path.isfile(pfile): continue
			openfunc = gzip.open if efile.endswith('.gz') else open
			with openfunc(pfile) as fpred, openfunc(efile) as fexpt:
				test.assertEqual(fpred.read().splitlines(), fexpt.read().splitlines())

def procOKIn(proc, msg, test):
	test.maxDiff = 1000
	predfile = proc.channel.get()
	if not isinstance(msg, list):
		msg = [msg]
	with open(predfile) as f:
		content = f.read()
		for m in msg:
			test.assertIn(m, content)

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

