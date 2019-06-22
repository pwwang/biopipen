from bioprocs.bin.arguments import commands

def params(opts):
	"""Show and query a parameter"""

def main():
	command, opts = commands._parse(arbi = True)
	if commands == 'params':
		params(opts)

if __name__ == "__main__":
	main()
