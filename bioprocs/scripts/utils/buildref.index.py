
from os import path, symlink, remove
from sys import stderr
if 'buildrefIndex' not in vars() or not callable (buildrefIndex):
	def buildrefIndex (jobid, infile1, outfiles1, cmd1, pdir, cmd2):
		
		def filesExist (files):
			for f in files:
				if not path.exists (f):
					return False
			return True
			
		if not isinstance (outfiles1, list):
			outfiles1 = [outfiles1]
			
		if filesExist (outfiles1): 
			return infile1
		
		infile2   = path.join (pdir, '0', 'output', path.basename(infile1))
		outfiles2 = [path.join(pdir, '0', 'output', path.basename(of1)) for of1 in outfiles1]
		
		if path.exists (infile2) and filesExist (outfiles2):
			stderr.write ("Index file found for '" + infile1 + "'.\n")
			return infile2
		
		donefile = outfiles2[0] + '.done'
		try:
			polling0 (jobid, cmd1, donefile)
			return infile1
		except SystemExit:
			if jobid == 0 and path.exists (donefile):
				remove (donefile)
			if jobid == 0 and path.exists (donefile + '.error'):
				remove (donefile)
			if not path.exists (infile2):
				symlink (infile1, infile2)
			polling0 (jobid, cmd2, donefile)
			return infile2