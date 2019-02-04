import helpers, testly
from collections import OrderedDict
from bioprocs.utils import shell

class TestFuncs(testly.TestCase):

	def dataProvider_testCmdargs(self):
		# 0: empty
		yield testly.Data(params = {}, ret = '')
		# 1: unexpected order
		yield testly.Data(params = dict(b = 1, a = 2), ret = '-a 2 -b 1')
		# 2: keeping order
		yield testly.Data(params = OrderedDict([('b', 1), ('a', 2)]), ret = '-b 1 -a 2')
		# 3: preceding positional arguments
		yield testly.Data(params = {'': 'x', 'b': 1, 'a': 2}, ret = 'x -a 2 -b 1')
		# 4: ending positional arguments
		yield testly.Data(params = {'_': 'x', 'b': 1, 'a': 2}, ret = '-a 2 -b 1 x')
		# 5: write stdout to a file
		yield testly.Data(params = {'_': 'x', 'b': 1, 'a': 2, '_stdout': 'file'}, ret = '-a 2 -b 1 x > file')
		# 6: append stdout to a file
		yield testly.Data(params = {'_': 'x', 'b': 1, 'a': 2, '__stdout': 'file'}, ret = '-a 2 -b 1 x >> file')
		# 7: dash and equal
		yield testly.Data(params = OrderedDict([('b', 1), ('ac', 2)]), ret = '-b 1 --ac=2')
		# 8: dash with '-'
		yield testly.Data(params = OrderedDict([('b', 1), ('ac', 2)]), dash = '-', ret = '-b 1 -ac=2')
		# 9: equal with ' '
		yield testly.Data(params = OrderedDict([('b', 1), ('ac', 2)]), equal = ' ', ret = '-b 1 --ac 2')
		# 10: bools, ignorefalse = True
		yield testly.Data(params = dict(a = False, b = True), ignorefalse = True, ret = '-b')
		# 11: bools, ignorefalse = False
		yield testly.Data(params = dict(a = False, b = True), ignorefalse = False, ret = '-a 0 -b')
		# 12: list
		yield testly.Data(params = dict(a = [1, 2, 3]), duplistkey = False, ret = '-a 1 2 3')
		# 13: list, duplistkey = True
		yield testly.Data(params = dict(a = [1, 2, 3]), duplistkey = True, ret = '-a 1 -a 2 -a 3')
		# 14: multiple same keys
		yield testly.Data(params = {'a': 1, 'a #1': 2, 'a #2': 3}, ret = '-a 1 -a 2 -a 3')
		# 15: quote
		yield testly.Data(params = dict(a = 'hello world', b = '"'), ret = '-a \'hello world\' -b \'"\'')

	def testCmdargs(self, params, ret, dash = 'auto', equal = 'auto', duplistkey = False, ignorefalse = True):
		self.assertEqual(shell.cmdargs(params, dash, equal, duplistkey, ignorefalse), ret)

	def dataProvider_testRuncmd(self):
		# 0: Normal
		yield testly.Data(cmd2run = 'ls', logs = ['RUNNING: ls', '- PID', '- STDERR', '- RETURNCODE: 0'])
		# 1: Stderr
		yield testly.Data(cmd2run = 'echo 1 1>&2; echo 2 1>&2', logs = ['RUNNING: echo', '- PID:', '- STDERR: 1', '- STDERR: 2', '- RETURNCODE: 0'])
		# 2: Exception
		yield testly.Data(cmd2run = 'nosuchcmd', raised = True)
		# 3: Another rc
		yield testly.Data(cmd2run = 'exit 127', rc = 127)

	def testRuncmd(self, cmd2run, raised = False, rc = 0, logs = None, **kwargs):
		if raised:
			with self.assertLogs('bioprocs', level = 'DEBUG'):
				self.assertRaises(shell.RuncmdException, shell.runcmd, cmd2run, raiseExc = True, **kwargs)
		else:
			with self.assertLogs('bioprocs', level = 'DEBUG') as cm:
				c = shell.runcmd(cmd2run, raiseExc = False, **kwargs)
			if logs:
				self.assertOrderedStrsInArray(logs, cm.output)
			self.assertEqual(c.rc, rc)

class TestShell(testly.TestCase):

	def dataProvider_testInit(self):
		yield None, {}
		yield object(), {'a': 1}

	def testInit(self, cmdobj, tools = None):
		sr = shell.ShellResult(cmdobj, tools)
		self.assertIs(sr.cmdobj, cmdobj)
		self.assertEqual(sr.tools, tools or {})
		self.assertIs(sr.done, False)



if __name__ == '__main__':
	testly.main(verbosity = 2)