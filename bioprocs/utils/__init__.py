
def alwaysList(l):
	ret = []
	if isinstance(l, (list, tuple)):
		for x in l:
			ret += alwaysList(x)
	else:
		ret = [x.strip() for x in l.split(',')]
	return ret