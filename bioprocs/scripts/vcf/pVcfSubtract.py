import gzip
from os import path, symlink, readlink, remove
from bioprocs.utils import runcmd, cmdargs
from subprocess import check_output
from bioprocs.utils.parallel import Parallel

infile1  = {{i.infile1 | quote}}
infile2  = {{i.infile2 | quote}}
outfile  = {{o.outfile | quote}}
rmany    = {{args.any | repr}}
tool     = {{args.tool | quote}}
header   = {{args.header | repr}}
bedtools = {{args.bedtools | quote}}
tabix    = {{args.tabix | quote}}
bychrom  = {{args.bychrom | repr}}
nthread  = int({{args.nthread | repr}})
outdir   = {{job.outdir | quote}}

if tool == 'bedtools':
	if header:
		# write the header to outfile
		openfunc = gzip.open if infile1.endswith('.gz') else open
		with openfunc(infile1) as fin, open(outfile, 'w') as fout:
			for line in fin:
				if not line.startswith('#'):
					break
				fout.write(line)

	# go directly
	if not bychrom:
		cmd = '{bedtools} subtract {params} >> {outfile}'
		params = {
			'a': infile1,
			'b': infile2
		}
		if rmany:
			params['A'] = True
		
		params = cmdargs(params, dash = '-', equal = ' ')
		runcmd(cmd.format(
			bedtools = bedtools,
			params   = params,
			outfile  = outfile
		))
	else:
		# check if infile1 is tabix indexed
		def tabindex(vcf, outdir):
			if not vcf.endswith('.gz'):
				gzfile = path.join(outdir, path.basename(vcf) + '.gz')
				bgzip_cmd = 'bgzip {!r} -c > {!r}'.format(
					vcf, 
					gzfile
				)
				runcmd(bgzip_cmd)
				runcmd('{} {!r}'.format(tabix, gzfile))
			else:
				gzfile = path.join(outdir, path.basename(vcf))
				# it is gzipped, try to find the index file (.tbi)
				idxfile = path.join(outdir, path.basename(vcf) + '.tbi')
				symlink(vcf, gzfile)
				while True:
					try:
						link = readlink(vcf)
						tbifile = link + '.tbi'
						if path.isfile(tbifile):
							symlink(tbifile, gzfile + '.tbi')
							break
						vcf = link
					except OSError:
						break
				
				if not path.isfile(idxfile):
					index_cmd = '{} {!r}'.format(tabix, gzfile)
					runcmd(index_cmd)

			return gzfile

		def runChrom(file1, file2, chrom):
			outfile1_list = list(file1[:-3].rpartition('.'))
			outfile1_list.insert(-2, '-' + chrom)
			outfile1 = ''.join(outfile1_list)
			outfile2 = list(file2[:-3].rpartition('.'))
			outfile2.insert(-2, '-' + chrom)
			outfile2 = ''.join(outfile2)
			outfile1_list.insert(-2, '.subtracted')
			outfile  = ''.join(outfile1_list)
			vfcmd    = '{} -h {!r} {!r} > {!r}'
			runcmd(vfcmd.format(tabix, file1, chrom, outfile1))
			runcmd(vfcmd.format(tabix, file2, chrom, outfile2))

			cmd = '{bedtools} subtract {params} > {outfile!r}'
			params = {
				'a': outfile1,
				'b': outfile2
			}
			if rmany:
				params['A'] = True
			
			params = cmdargs(params, dash = '-', equal = ' ')
			runcmd(cmd.format(
				bedtools = bedtools,
				params   = params,
				outfile  = outfile
			))
			remove(outfile1)
			remove(outfile2)
			return outfile
		
		infile1 = tabindex(infile1, outdir)
		infile2 = tabindex(infile2, outdir)
		chroms  = [chr.strip() for chr in check_output([tabix, '-l', infile1]).splitlines()]

		if nthread > 1:
			p = Parallel(nthread, raiseExc = True)
			outfiles = p.run(runChrom, [(infile1, infile2, chrom) for chrom in chroms])
			# make sure it's in the right order
			outfiles = sorted(outfiles, key = lambda x: chroms.index(x.split('.')[-2]))
		else:
			outfiles = []
			for chrom in chroms:
				outfiles.append(runChrom(infile1, infile2, chrom))

		with open(outfile, 'a+') as fout:
			for of in outfiles:
				with open(of) as f:
					fout.write(f.read())
		
elif tool == 'pyvcf':
	import vcf
	v1  = vcf.Reader(filename = infile1)
	v2  = vcf.Reader(filename = infile2)
	out = vcf.Writer(open(outfile, 'w'), v1)
	r1  = r2 = None
	while True:
		try:
			r1 = r1 or v1.next()
			r2 = r2 or v2.next()
			if r1.CHROM < r2.CHROM:
				out.write_record(r1)
				r1 = None
				continue
			if r1.CHROM > r2.CHROM:
				r2 = None
				continue
			if r1.POS < r2.POS:
				out.write_record(r1)
				r1 = None
				continue
			if r1.POS > r2.POS:
				r2 = None
				continue
			r1 = None
		except StopIteration:
			break
		except:
			r1 = r2 = None
			continue
else:
	import gzip
	with open(infile1) if not infile1.endswith('.gz') else gzip.open(infile1) as f1, \
		 open(infile2) if not infile2.endswith('.gz') else gzip.open(infile2) as f2, \
		 open(outfile, 'w') as fout:
		line1 = line2 = None
		for line in f1:
			if line.startswith('#'):
				fout.write(line)
			else:
				line1 = line.split('\t')
				break
		
		while True:
			try:
				line1 = line1 or f1.next().split('\t')
				line2 = line2 or f2.next().split('\t')
				if line1[0] < line2[0]:
					fout.write('\t'.join(line1))
					line1 = None
					continue
				if line1[0] > line2[0]:
					line2 = None
					continue
				if int(line1[1]) < int(line2[1]):
					fout.write('\t'.join(line1))
					line1 = None
					continue
				if int(line1[1]) > int(line2[1]):
					line2 = None
					continue
				line1 = None
			except StopIteration:
				break
			except:
				line1 = line2 = None
				continue
			


