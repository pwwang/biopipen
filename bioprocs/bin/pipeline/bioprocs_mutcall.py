"""Call mutatoins from 2nd-gen sequencing data"""

import sys
from pyppl import Box
from bioprocs import params


def main(args = None, prog = None):
	prog     = prog or sys.argv[0]
	args     = args or sys.argv[1:]
	params._assembler.prog = prog
	params._parse(args, dict_wrapper = Box)
	print (params)

if __name__ == "__main__":
	main()