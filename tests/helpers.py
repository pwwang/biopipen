import testly

def assertOrderedStrsInArray(self, strs, array, msg = None):
	if not self.maxDiff is None:
		self.maxDiff = max(self.maxDiff or 5000, 5000)
	
	self.assertIsInstance(strs,  list, 'First argument is not a list.')
	self.assertIsInstance(array, list, 'Second argument is not a list.')

	for s in strs:
		while array:
			a = array.pop(0)
			if s in a:
				array.append(None)
				break
			continue
		if array: # found
			continue
		standardMsg = '%s not in %s' % (testly.util.safe_repr(s, True), testly.util.safe_repr(array, True))
		self.fail(self._formatMessage(msg, standardMsg))

testly.TestCase.assertOrderedStrsInArray = assertOrderedStrsInArray