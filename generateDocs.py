"""
@description:
	Generate documentation for bioprocs, if doctoc is installed, try to generate TOC.
@sources:
	bioprocs/*.py (__init__.py not included)
@destination:
	docs/DOCS.md
@fields:
	- name: required, the name of the process (must follow \""" or ''')
	- description: description of the process
	- input: input of the process
	- output: output of the process
	- args: args could be set when call the process
	- requires: tools required by this process
"""
import re, glob, os, bioprocs

srcdir  = os.path.join (  os.path.dirname (os.path.realpath(__file__)), 'bioprocs')
dstfile = os.path.join (  os.path.dirname (os.path.realpath(__file__)), 'docs', 'DOCS.md')
infiles = [sfile for sfile in glob.glob(os.path.join(srcdir, '*.py')) if os.path.basename (sfile) != "__init__.py" ]

fout    = open (dstfile, 'w')

def header():
	fout.write ("""
# Documentation for bioprocs v%s
A set of procs for bioinformatics using [pyppl](https://github.com/pwwang/pyppl)
""" % bioprocs.VERSION)

def each (infile):
	name = os.path.basename(infile)[:-3].upper()
	fout.write ("\n## %s\n" % name)
	blocks = re.findall (r'((\"\"\"|\'\'\')\s*\n@name:[\s\S]+?\2)', open(infile).read())
	for block in blocks:
		block = block[0][3:-3]
		lines = block.split("\n")
		fieldname = ''
		content   = ''
		lines.append ('@')
		for i in range(len(lines)):
			line = lines[i]
			if line.startswith ("@"):
				if fieldname == 'name':
					fout.write ("\n### %s\n" % content)
				elif fieldname:
					fout.write ("#### %s\n" % fieldname)
					fout.write ("%s\n" % content)

				fields    = line[1:].split(':')
				fieldname = fields[0].strip()

				if len(fields) > 1:
					c = ''.join(fields[1:]).strip()
					if not c:
						content = ''
					else:
						if fieldname == 'name':
							content = ''.join(fields[1:])
						else:
							content = '-' + ''.join (fields[1:]) + '\n'
			else:
				if fieldname == 'name':
					content += ' ' + line.lstrip("\t")
				else:
					l = line.lstrip(" \t")
					if not l: continue
					content += '- ' + l + '\n'


def doctoc ():
	pass



def main ():
	header()
	for infile in infiles:
		each (infile)

	fout.close()
	doctoc()


if __name__ == "__main__":
	main ()


