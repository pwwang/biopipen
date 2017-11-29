"""
@description:
	Generate documentation for bioprocs, if doctoc is installed, try to generate TOC.
@sources:
	bioprocs/*.py (__init__.py not included)
@destination:
	docs/*.md
@fields:
	- name: required, the name of the process (must follow \""" or ''')
	- description: description of the process
	- input: input of the process
	- output: output of the process
	- args: args could be set when call the process
	- requires: tools required by this process
"""
import re, bioprocs
from sys import stdout
from glob import glob
from os import path, symlink, devnull

srcdir  = path.join (path.dirname (path.realpath(__file__)), 'bioprocs')
dstdir  = path.join (path.dirname (path.realpath(__file__)), 'docs')
infiles = [f for f in glob(path.join(srcdir, '*.py')) if path.basename(f) != "__init__.py" ]

def summary():
	summaryfile = path.join(dstdir, 'SUMMARY.md')
	readmefile  = path.join(dstdir, 'README.md')
	if not path.exists(readmefile):
		symlink(path.join('../', 'README.md'), readmefile)

	with open(summaryfile, 'w') as fout:
		fout.write('# Summary\n\n')
		fout.write('* [Introduction](README.md)\n')
		for infile in infiles:
			title = path.basename(infile)[:-3]
			fout.write('* [%s](%s.md)\n' % (title, title))

def eachfile(infile):
	name   = path.basename(infile)[:-3]
	mdfile = path.join(dstdir, path.basename(infile)[:-3] + '.md')
	with open(infile) as f, open(mdfile, 'w') as fout:
		fout.write('# %s\n' % name)
		fout.write('<!-- toc -->\n')

		fout.write('{% raw %}\n')
		content = f.read()
		# docs for module?
		moddocs = re.match(r'^(\"\"\"|\'\'\')\s*\n(@\w+:[\s\S]+?)\1', content)
		if moddocs:
			modsec = moddocs.group(1)
			fout.write(fmtSection(modsec, True))
			content = content[len(modsec):]
		
		blocks = re.findall (r'(\"\"\"|\'\'\')\s*\n(@name:[\s\S]+?)\1', content)
		for block in blocks:
			block   = block[1].splitlines()
			section = []
			for line in block:
				if line.startswith('@'):
					if section:
						fout.write(fmtSection(section))
					section = [line]
				else:
					section.append(line)
		fout.write('{% endraw %}\n')


def fmtSection(lines, mod = False):
	title = lines.pop(0).strip()
	basepunk = '#' * (1 if mod else 2)
	if title.startswith('@name'):
		return '\n%s %s\n' % (basepunk, lines[0].strip())
	
	ret = '\n%s# %s\n' % (basepunk, title.strip('@:'))
	for line in lines:
		secname =  re.match('^\t(`.+?`\s*:)', line)
		if secname:
			secname = secname.group(1)
			ret += '%s## %s\n' % (basepunk, secname)
			ret += line[len(secname)+1:].lstrip() + '  \n'
		else:
			ret += line[1:] + '\n'
	return ret

if __name__ == "__main__":
	summary()
	for infile in infiles:
		stdout.write('- Compiling %s\n' % infile)
		eachfile (infile)
	stdout.write('Done!\n')

