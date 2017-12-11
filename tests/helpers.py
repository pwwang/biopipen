import inspect, gzip
from subprocess import Popen, PIPE
from pyppl import PyPPL, Proc
from os import path, listdir

filedir = None
config  = {'log': {'level': 'basic', 'lvldiff': ["+P.ARGS", "-DEBUG"]}}

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

def fileOK(predfile, exptfile, test, order=True, comment=None):
	openfunc = gzip.open if exptfile.endswith('.gz') else open
	with openfunc(predfile) as fpred, openfunc(exptfile) as fexpt:
		predlines = [line for line in fpred.read().splitlines() if not comment or not line.startswith(str(comment))]
		exptlines = [line for line in fexpt.read().splitlines() if not comment or not line.startswith(str(comment))]
		if order:
			test.assertEqual(predlines, exptlines)
		else:
			test.assertEqual(sorted(predlines), sorted(exptlines))

def fileOKIn(exptfile, msg, test):
	openfunc = gzip.open if exptfile.endswith('.gz') else open
	if not isinstance(msg, list):
		msg = [msg]
	with openfunc(exptfile) as f:
		content = f.read()
		for m in msg:
			test.assertIn(m, content)

def procOK(proc, name, test, order = True, comment = None):
	predfile = proc.channel.get()
	exptfile = getfile(name, False)
	if path.isfile(predfile):
		fileOK(predfile, exptfile, test, order, comment)
	else:
		for item in listdir(predfile):
			pfile = path.join(predfile, item)
			efile = path.join(exptfile, item)
			if not path.isfile(pfile): continue
			fileOK(pfile, efile, test, order, comment)
			
def procOKIn(proc, msg, test):
	fileOKIn(proc.channel.get(), msg, test)

def cmdOK(cmd, test, inout = None, inerr = None, testrc = True):
	p = Popen(cmd, stdout = PIPE, stderr = PIPE)
	out, err = p.communicate()
	if inout:
		if not isinstance(inout, list):
			inout = [inout]
		for ino in inout:
			test.assertIn(ino, out)
	if inerr:
		if not isinstance(inerr, list):
			inerr = [inerr]
		for ine in inerr:
			test.assertIn(ine, err)
	if testrc:
		rc = p.returncode
		test.assertEqual(rc, 0)

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

