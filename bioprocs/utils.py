from pyppl import proc, doct

mem = doct({
	'toJava': {},
	'toM': {},
})
mem.toJava.python = """
if 'memtoJava' not in vars() or not callable (memtoJava):
	def memtoJava (mem):
		if mem.endswith('G') or mem.endswith('g'):
			num = int (mem[:-1])
			if num >= 8:
				xms = "-Xms%sG" % str(int (num/8))
			else:
				xms = "-Xms%sM" % str (num * 128)
			return xms + " -Xmx%sG" % str(num)
		elif mem.endswith('M') or mem.endswith('m'):
			num = int (mem[:-1])
			xms = "-Xms%sM" % str (num/8)
			return xms + " -Xmx%sM" % str(num)
		else:
			import sys
			sys.stderr.write ('Unknown unit of mem, expect G/g/M/m')
			sys.exit (1)
"""

mem.toJava.r = """
if (!exists('memtoJava')) {
	memtoJava = function (mem) {
		if (grepl('G$', mem) || grepl('g$', mem)) {
			num = as.numeric (substr(mem, 1, nchar(mem)-1))
			if (num >= 8) {
				xms = paste('-Xms', num/8, 'G', sep='')
			} else {
				xms = paste('-Xms', num*128, 'M', sep='')
			}
			return (paste(xms, ' -Xmx', num, 'G', sep=''))
		} else if (grepl('M$', mem) || grepl('m$', mem)) {
			num = as.numeric (substr(mem, 1, nchar(mem)-1))
			return (paste('-Xms', num/8, 'M', ' -Xmx', num, 'M', sep=''))
		} else {
			stop ('Unknown unit of mem, expect G/g/M/m')
		}
		
	}
}
"""

mem.toM.python = """
if 'memtoM' not in vars() or not callable(memtoM):
	def memtoM (mem):
		if mem.endswith('G') or mem.endswith('g'):
			num = int (mem[:-1])
			return str (num * 1024)
		elif mem.endswith('M') or mem.endswith('m'):
			return mem[:-1]
		else:
			import sys
			sys.stderr.write ('Unknown unit of mem, expect G/g/M/m')
			sys.exit (1)
"""

runcmd = doct()
runcmd.python = """
from subprocess import Popen
from sys import stderr, exit as SysExit
if 'runcmd' not in vars() or not callable(runcmd):
	def runcmd (cmd, exit = True):
		stderr.write ("%s\\n" % ('-'*80))
		stderr.write ("Running command:\\n")
		stderr.write (cmd + "\\n")
		rc = Popen (cmd, shell=True).wait()
		if rc == 0:
			stderr.write ("Return code: 0\\n")
			stderr.write ("%s\\n" % ('-'*80))
		else:
			stderr.write ("Return code: %s\\n" % str(rc))
			stderr.write ("%s\\n" % ('-'*80))
			if exit:
				SysExit(1)
		return rc == 0
"""

runcmd.r = """
if (!exists('runcmd')) {
	runcmd = function (cmd, exit = TRUE, intern=FALSE) {
		write (paste(rep('-', 80), collapse='') , stderr())
		write ("Running command: ", stderr())
		write (cmd, stderr())
		
		if (!intern) {
			rc = system (cmd)
			if (rc == 0) {
				write ("Return code: 0", stderr())
				write (paste(rep('-', 80), collapse='') , stderr())
			} else {
				write (paste("Return code: ", rc, sep=''), stderr())
				write (paste(rep('-', 80), collapse='') , stderr())
				if (exit) {
					stop ('stop.')
				}
			}
			return (rc)
		} else {
			return (system(cmd, intern=T))
		}
	}
}
"""

buildArgsFastaFai = doct()
buildArgsFastaFai.bash   = """
ref="{{args.ref}}"
fai="$ref.fai"
if [[ ! -e "$fai" ]]; then
	echo "Generating fai index for $ref ..." 1>&2
	{{args.samtools}} faidx "$ref"
fi
"""
buildArgsFastaDict = doct()
buildArgsFastaDict.bash   = """
ref="{{args.ref}}"
dict="${ref%.*}.dict"
if [[ ! -e "$dict" ]]; then
	echo "Generating dict index for $ref ..." 1>&2
	{{args.picard}} CreateSequenceDictionary R="$ref" O="$dict"
fi
"""
checkArgsRef = doct()
checkArgsRef.bash = """
if [[ ! -e "{{args.ref}}" ]]; then
	echo "Reference file is not specfied or not exists." 1>&2
	exit 1
fi
"""
checkArgsRefgene = doct()
checkArgsRefgene.bash = """
if [[ ! -e "{{args.refgene}}" ]]; then
	echo "Reference gene file is not specfied or not exists." 1>&2
	exit 1
fi
"""

polling0 = doct()
polling0.python = """
%s

from os import path
from time import sleep
from sys import stderr, exit

if 'polling0' not in vars() or not callable(polling0):
	def polling0 (jobid, cmd, flagfile, t = 10):
		errorfile = flagfile + '.error'
		
		if jobid == 0:
			try:
				runcmd (cmd)
				open (flagfile, 'w').close()
			except SystemExit:
				open (errorfile, 'w').close()
		else:
			while True:
				if path.exists(errorfile):
					stderr.write('Error happened for job #0!\\n')
					exit (1)
				if path.exists(flagfile):
					return
				stderr.write('Waiting for job #0 ...\\n')
				sleep(t)
				
""" % (runcmd.python)

# requires all jobs run simultaneously
pollingAll = doct()
pollingAll.python = """
%s

from os import path
from time import sleep
from sys import stderr, exit

if 'pollingAll' not in vars() or not callable(pollingAll):
	def pollingAll (workdir, length, jobid, cmd, flagfname, t = 10):
		errorfname = flagfname + '.error'
		flagfile   = path.join(workdir, str(jobid), 'output', flagfname)
		errorfile  = path.join(workdir, str(jobid), 'output', errorfname)
		
		try:
			runcmd (cmd)
			open (flagfile, 'w').close()
		except SystemExit:
			open (errorfile, 'w').close()
			
		wait = True
		while wait:
			wait = False
			for i in range(length):
				efile = path.join(workdir, str(i), 'output', errorfname)
				ffile = path.join(workdir, str(i), 'output', flagfname)
				if path.exists (efile):
					stderr.write ('Job '+ str(i) +' failed, I am also exiting ...')
					exit (1)
				if not path.exists (ffile):
					wait = True
					break
			if wait:
				stderr.write('Waiting till jobs done ...\\n')
				sleep(t)
			else:
				break
		
				
""" % (runcmd.python)

pollingAll.r = """
%s

if (!exists('pollingAll')) {
	pollingAll = function (workdir, length, jobid, cmd, flagfname, t = 10) {
		errorfname = paste(flagfname, '.error', sep='')
		flagfile   = file.path (workdir, jobid, 'output', flagfname)
		errorfile  = file.path (workdir, jobid, 'output', errorfname)
		
		tryCatch ({
			runcmd (cmd)
			file.create (flagfile)
		}, error = function(cond) {
			file.create (errorfile)
		})
		
		wait = TRUE
		while (wait) {
			wait = FALSE
			for (i in 0:(length-1)) {
				efile = file.path (workdir, i, 'output', errorfname)
				ffile = file.path (workdir, i, 'output', flagfname)
				
				if (file.exists(efile)) {
					stop (paste('Job', i, 'failed, I am also exiting ...'), stderr())
				}
				if (!file.exists(ffile)) {
					write (paste('file not exists:', ffile), stderr())
					wait = TRUE
					break
				}
			}
			if (wait) {
				write('Waiting till all jobs done ...\\n', stderr())
				Sys.sleep(t)
			} else {
				break
			}
		}
		
	}
}
""" % (runcmd.r)

buildArgIndex = doct ()
buildArgIndex.python = """
%s

from os import path, symlink, remove
from sys import stderr
if 'buildArgIndex' not in vars() or not callable (buildArgIndex):
	def buildArgIndex (jobid, infile1, outfiles1, cmd1, pdir, cmd2):
		
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
			stderr.write ("Index file found for '" + infile1 + "'.\\n")
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
		
""" % (polling0.python)

cbindFill = doct()
cbindFill.r = """
if (!exists('cbindFill')) {
	cbindFill = function (x1, x2) {
		y = merge(x1, x2, by='row.names', all=T, sort=F)
		rownames(y) = y[, "Row.names"]
		y = y[, -1, drop=F]
		cnames      = c(colnames(x1), colnames(x2))
		if (!is.null(cnames)) {
			colnames(y) = cnames
		}
		return (y)
	}
}
"""

plot = doct({
	'hist':    {},
	'heatmap': {},
	'symbols': {},
	'text':    {},
	'venn':    {},
	'upset':   {},
})
plot.hist.r = """
require ('ggplot2')
if (!exists('plotHist')) {
	plotHist = function (m, filename, xlab = " ", ylab="Frequency", main=" ") {
		png (filename, res=300, width=2000, height=2000)
		print (
			ggplot(stack(as.data.frame(t(m))))  +
				geom_histogram(aes(values), stat='count') +
				ggtitle(main) +
				scale_y_continuous(name = ylab) +
				scale_x_continuous(name = xlab) +
				theme(
					plot.title = element_text(hjust = 0.5),
					axis.text.x = element_text(angle = 60, hjust = 1)
				)
			
		)
		dev.off()
	}
}
"""

plot.heatmap.r = """
require('pheatmap')
if (!exists('plotHeatmap')) {
	plotHeatmap = function (mat, params = list()) {
		do.call(pheatmap, c(list(mat=mat), params))
	}
}
"""

plot.symbols.r = """
if (!exists('plotSymbols')) {
	plotSymbols = function (x, y, params) {
		do.call(symbols, c(list(x=x, y=y), params))
	}
}
"""

plot.text.r = """
if (!exists('plotText')) {
	plotText = function (x, y, params) {
		do.call(text, c(list(x=x, y=y), params))
	}
}
"""

plot.venn.r = """
if (!exists('plotVenn')) {
	plotVenn = function(mat, filename, params) {
		library(VennDiagram)
		nc   = ncol(mat)
		cats = colnames(mat)
		csum = colSums(mat)
		rsum = rowSums(mat)
		png (filename, res=300, width=2000, height=2000)
		if (nc == 1) {
			v = do.call(draw.single.venn, c(list(area = csum[1,1], category = cats[1]), params))
		} else if (nc == 2) {
			cross.area = length(which(rsum == nc))
			v = do.call(draw.pairwise.venn, c(list(area1 = csum[1], area2 = csum[2], cross.area = cross.area, category = cats), params))
		} else if (nc == 3) {
			n12  = length(which(rowSums(mat[, c(1, 2)]) == 2))
			n23  = length(which(rowSums(mat[, c(2, 3)]) == 2))
			n13  = length(which(rowSums(mat[, c(1, 3)]) == 2))
			n123 = length(which(rsum == nc))
			v = do.call(draw.triple.venn, c(list(area1 = csum[1], area2 = csum[2], area3 = csum[3], n12 = n12, n23 = n23, n13 = n13, n123 = n123, category = cats), params))
		} else if (nc == 4) {
			n12  = length(which(rowSums(mat[, c(1, 2)]) == 2))
			n13  = length(which(rowSums(mat[, c(1, 3)]) == 2))
			n14  = length(which(rowSums(mat[, c(1, 4)]) == 2))
			n23  = length(which(rowSums(mat[, c(2, 3)]) == 2))
			n24  = length(which(rowSums(mat[, c(2, 4)]) == 2))
			n34  = length(which(rowSums(mat[, c(3, 4)]) == 2))
			n123 = length(which(rowSums(mat[, c(1,2,3)]) == 3))
			n124 = length(which(rowSums(mat[, c(1,2,4)]) == 3))
			n134 = length(which(rowSums(mat[, c(1,3,4)]) == 3))
			n234 = length(which(rowSums(mat[, c(2,3,4)]) == 3))
			n1234 = length(which(rsum == nc))
			v = do.call(draw.quad.venn, c(list(
				area1 = csum[1], 
				area2 = csum[2], 
				area3 = csum[3], 
				area4 = csum[4], 
				n12   = n12, 
				n13   = n13, 
				n14   = n14, 
				n23   = n23, 
				n24   = n24, 
				n34   = n34, 
				n123  = n123, 
				n124  = n124, 
				n134  = n134, 
				n234  = n234, 
				n1234 = n1234, 
				category = cats), params))
		} else if (nc == 5) {
			n12    = length(which(rowSums(mat[, c(1, 2)]) == 2))
			n13    = length(which(rowSums(mat[, c(1, 3)]) == 2))
			n14    = length(which(rowSums(mat[, c(1, 4)]) == 2))
			n15    = length(which(rowSums(mat[, c(1, 5)]) == 2))
			n23    = length(which(rowSums(mat[, c(2, 3)]) == 2))
			n24    = length(which(rowSums(mat[, c(2, 4)]) == 2))
			n25    = length(which(rowSums(mat[, c(2, 5)]) == 2))
			n34    = length(which(rowSums(mat[, c(3, 4)]) == 2))
			n35    = length(which(rowSums(mat[, c(3, 5)]) == 2))
			n45    = length(which(rowSums(mat[, c(4, 5)]) == 2))
			n123   = length(which(rowSums(mat[, c(1, 2, 3)]) == 3))
			n124   = length(which(rowSums(mat[, c(1, 2, 4)]) == 3))
			n125   = length(which(rowSums(mat[, c(1, 2, 5)]) == 3))
			n134   = length(which(rowSums(mat[, c(1, 3, 4)]) == 3))
			n135   = length(which(rowSums(mat[, c(1, 3, 5)]) == 3))
			n145   = length(which(rowSums(mat[, c(1, 4, 5)]) == 3))
			n234   = length(which(rowSums(mat[, c(2, 3, 4)]) == 3))
			n235   = length(which(rowSums(mat[, c(2, 3, 5)]) == 3))
			n245   = length(which(rowSums(mat[, c(2, 4, 5)]) == 3))
			n345   = length(which(rowSums(mat[, c(3, 4, 5)]) == 3))
			n1234  = length(which(rowSums(mat[, c(1, 2, 3, 4)]) == 4))
			n1235  = length(which(rowSums(mat[, c(1, 2, 3, 5)]) == 4))
			n1245  = length(which(rowSums(mat[, c(1, 2, 4, 5)]) == 4))
			n1345  = length(which(rowSums(mat[, c(1, 3, 4, 5)]) == 4))
			n2345  = length(which(rowSums(mat[, c(2, 3, 4, 5)]) == 4))
			n12345 = length(which(rsum == nc))
			v = do.call(draw.quintuple.venn, c(list(
				area1  = csum[1], 
				area2  = csum[2], 
				area3  = csum[3], 
				area4  = csum[4], 
				area5  = csum[5], 
				n12    = n12, 
				n13    = n13, 
				n14    = n14, 
				n15    = n15, 
				n23    = n23, 
				n24    = n24, 
				n25    = n25, 
				n34    = n34, 
				n35    = n35, 
				n45    = n45, 
				n123   = n123, 
				n124   = n124, 
				n125   = n125, 
				n134   = n134, 
				n135   = n135, 
				n145   = n145, 
				n234   = n234, 
				n235   = n235, 
				n245   = n245, 
				n345   = n345, 
				n1234  = n1234, 
				n1235  = n1235, 
				n1245  = n1245, 
				n1345  = n1345, 
				n2345  = n2345, 
				n12345 = n12345, 
				category = cats), params))
		} else {
			stop('Too many columns to draw a venn diagram, please use UpSetR.')
		}
		print (v)
		dev.off()
	}
}
"""

plot.upset.r = """
if (!exists('plotUpset')) {
	plotUpset = function(mat, filename, params) {
		library(UpSetR)
		png (filename, res=300, width=2000, height=2000)
		if (! "nintersects" %in% names(params)) {
			params$nintersects = NA
		}
		v = do.call (upset, c(list(data = mat, sets = colnames(mat)), params))
		print (v)
		dev.off()
	}
}
"""

txt = doct({
	'filter': {},
})
txt.filter.python = """
if 'txtFilter' not in vars() or not callable (txtFilter):
	def txtFilter(infile, outfile, cols, rfilter, header = True, skip = 0, delimit = "\\t"):
		if not isinstance(cols, list):
			cols = map(lambda x: x.strip(), cols.split(','))
			cols = [int(c) if c.isdigit() else c for c in cols]
		
		if not cols and not rfilter and not skip:
			import os
			os.rename(infile, outfile)
		
		import csv
		with open (infile, 'r') as f, open(outfile, 'w') as fout:
			fcsv = csv.reader(f, delimiter = delimit)
			for x in range(skip):
				next(fcsv)
				
			i = 0
			for parts in fcsv:
				if i == 0:
					i += 1
					if header:
						cols = [c if isinstance(c, int) else parts.index(c) for c in cols]
				else:
					if rfilter and not rfilter(parts):
						continue
				fout.write("\\t".join([h for i,h in enumerate(parts) if i in cols]) + "\\n")
		
"""

download = doct({
	'curl': {},
})
download.curl.python = """
%s

if 'downloadCurl' not in vars() or not callable(downloadCurl):
	def downloadCurl(url, outfile, username = '', password = '', curl = 'curl'):
		"curl -Lfv -o filename.zip -u EMAIL:PASSWORD https://www.drugbank.ca/releases/5-0-7/downloads/target-approved-uniprot-links"
		cmd  = curl + ' -L '
		if username and password:
			cmd += '-u {}:{}'.format(username, password)
		elif username:
			cmd += '-u {}'.format(username)
		cmd += ' -o "{}" "{}"'.format(outfile, url)
		runcmd (cmd)
		
""" % (runcmd.python)