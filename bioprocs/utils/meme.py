"""
Read or write MEME motif file
"""
import re
import math
from collections import OrderedDict

class MemeRecord(object):

	def __init__(self, 
		name, 
		matrix,
		altname = '',
		mtype   = 'letter-probability',
		alength = None,
		w       = None,
		nsites  = 20,
		E       = 0,
		URL     = None
	):
		self.name    = name
		self.matrix  = matrix
		self.altname = altname
		self.mtype   = mtype
		self.alength = alength or len(matrix[0])
		self.w       = w or len(matrix)
		self.nsites  = nsites
		self.E       = E
		self.URL     = URL

	def __str__(self):
		return """
MOTIF {name}{altname}
{mtype} matrix: alength= {alength} w= {w} nsites= {nsites} E= {E}
{matrix}
{URL}
""".format(
		name    = self.name,
		altname = " " + self.altname if self.altname else "",
		mtype   = self.mtype,
		alength = self.alength,
		w       = self.w,
		nsites  = self.nsites,
		E       = self.E,
		matrix  = "\n".join(" ".join(str(r) for r in row) for row in self.matrix),
		URL     = "URL {}".format(self.URL) if self.URL else ""
	)

	def pwm2logodds(self):
		assert self.mtype == 'letter-probability'
		matrix = [
			tuple(math.exp(p)/(1.0 + math.exp(p)) for p in row)
			for row in self.matrix
		]
		return MemeRecord(
			name    = self.name,
			matrix  = matrix,
			altname = self.altname,
			mtype   = 'log-odds',
			alength = self.alength,
			w       = self.w,
			nsites  = self.nsites,
			E       = self.E,
			URL     = self.URL
		)

	def pwm2prob(self):
		assert self.mtype == 'log-odds'
		matrix = [
			tuple(math.log(p/(1.0-p)) for p in row)
			for row in self.matrix
		]
		return MemeRecord(
			name    = self.name,
			matrix  = matrix,
			altname = self.altname,
			mtype   = 'letter-probability',
			alength = self.alength,
			w       = self.w,
			nsites  = self.nsites,
			E       = self.E,
			URL     = self.URL
		)

class MemeReader(object):

	def __init__(self, memefile):
		self.meta     = {}
		alphabet_flag = False
		bgfreqs_flag  = False
		self.file = open(memefile)
		self.tell = 0
		while True:
			self.tell = self.file.tell()
			line = self.file.readline()
			if not line:
				raise ValueError('Not a valid MEME motif file.')
			if line.startswith('MEME version'):
				self.meta['Version'] = line[12:].strip()
			elif line.startswith('ALPHABET='):
				self.meta['Alphabet'] = line[9:].strip()
			elif line.startswith('ALPHABET'):
				self.meta['Alphabet'] = line[8:].strip()
				alphabet_flag = True
			elif line.startswith('END ALPHABET'):
				alphabet_flag = False
			elif alphabet_flag:
				self.meta['Alphabet'] += '\n' + line.strip()
			elif line.startswith('strands:'):
				self.meta['Strands'] = line[8:].strip()
			elif line.startswith('Background letter frequencies'):
				bgfreqs_flag = True
				source = line[30:].strip()
				if source.startswith('(from '):
					source = source[6:-2]
				else:
					source = ''
				self.meta['bgfreqs'] = {'from': source, 'freqs': OrderedDict()}
			elif bgfreqs_flag:
				bgfreqs_flag = False
				parts = line.strip().split()
				self.meta['bgfreqs']['freqs'] = OrderedDict(tuple([parts[2*i], float(parts[2*i+1])] for i  in range(int(len(parts)/2))))
			elif line.startswith('MOTIF'):
				self.file.seek(self.tell)
				break

	def next(self):
		name    = None
		altname = ''
		url     = None
		mtype   = ''
		matrix  = []
		attrs   = {}
		while True:
			tell = self.file.tell()
			line = self.file.readline()
			if not line:
				raise StopIteration()
			if line.startswith('MOTIF'):
				if name:
					self.file.seek(tell)
					break
				parts = line[5:].strip().split()
				name  = parts.pop(0)
				if parts:
					altname = parts[0]
			elif line.startswith('URL'):
				url = line[3:].strip()
			elif 'matrix:' in line:
				matrix = [] # in case there are multiple matrices
				mtype, attrs = line.strip().split('matrix:')
				mtype = mtype.strip()
				attrs = re.split(r'(?:\s*=\s*|\s+)', attrs.strip())
				attrs = {attrs[2*i]:attrs[2*i+1] for i in range(int(len(attrs)/2))}
			else:
				line = line.strip()
				if not line:
					continue
				matrix.append(tuple(float(v) for v in line.split()))
		
		return MemeRecord(
			name, 
			matrix,
			altname = altname,
			mtype   = mtype,
			URL     = url,
			**attrs
		)
	
	def __next__(self):
		return self.next()

	def rewind(self):
		self.file.seek(self.tell)

	def __iter__(self):
		return self

	def __del__(self):
		self.close()

	def close(self):
		if self.file:
			self.file.close()

class MemeWriter(object):

	def __init__(self, outfile, meta = None):
		self.meta = meta or {}
		self.file = open(outfile, 'w')

	def writeMeta(self):
		self.file.write("MEME version {}\n\n".format(self.meta.get('Version', 4)))
		alphabet = self.meta.get('Alphabet', 'ACGT')
		if '\n' in alphabet:
			self.file.write("ALPHABET {}\nEND ALPHABET\n\n".format(alphabet))
		else:
			self.file.write("ALPHABET= {}\n\n".format(alphabet))
		strands = self.meta.get('Strands', '+ -')
		self.file.write("strands: {}\n\n".format(strands))
		bgfreqs = self.meta.get("bgfreqs", {})
		if "from" in bgfreqs:
			self.file.write("Background letter frequencies (from {}):\n".format(bgfreqs['from']))
		if "freqs" in bgfreqs:
			self.file.write(" ".join('{} {}'.format(k, v) for k, v in bgfreqs['freqs'].items()) + "\n\n")
	
	def write(self, mrec):
		self.file.write(str(mrec))

	def __del__(self):
		self.close()

	def close(self):
		if self.file:
			self.file.close()
