from shutil import move, rmtree
from os import makedirs, path, symlink, remove
from sys import stdout, stderr
from collections import OrderedDict

{{ runcmd }}
{{ mem2 }}
{{ params2CmdArgs }}

infile    = {{ in.infile | quote }}
outfile   = {{ out.outfile | quote }}
informat  = {{ args.infmt | quote }}
informat  = informat if informat else {{ in.infile | ext | [1:] | quote}} 
tmpdir    = path.join ("{{args.tmpdir}}", "{{proc.id}}.{{in.infile | fn}}.{{job.index}}")
doSort    = {{ args.sort }}
doIndex   = {{ args.index }}
doMarkdup = {{ args.markdup }}
doRmdup   = {{ args.rmdup }}
params    = OrderedDict({{ args.params }})
if doRmdup:
	doMarkdup = True
if not path.exists (tmpdir):
	makedirs (tmpdir)
try:
	############### biobambam
	{% if args.tool | lambda x: x == 'biobambam' %}
	mem = mem2({{ args.mem | quote }}, 'M')
	if doIndex:
		params['index'] = 1
		params['indexfilename'] = {{out.idxfile | quote}}
	params['I']              = infile
	params['O']              = outfile
	params['SO']             = {{args.sortby | quote}}
	params['blockme']        = mem
	params['tmpfile']        = path.join(tmpdir, 'tmp.')
	params['inputformat']    = informat
	params['outfmt']         = 'bam'
	params['inputthreads']   = {{args.nthread}}
	params['outputthreads']  = {{args.nthread}}
	params['markduplicates'] = int(doMarkdup)
	params['rmdup']          = int(doRmdup)
	

	cmd = '{{args.biobambam_bamsort}} %s' % params2CmdArgs(params, dash = '', equal = '=', noq = ['index', 'inputthreads', 'outputthreads', 'markduplicates', 'rmdup'])
	runcmd (cmd)

	############### sambamba
	{% elif args.tool | lambda x: x == 'sambamba' %}
	if not (doSort or doIndex or doMarkdup or doRmdup):
		cmd = '{{args.sambamba}} view -S -f bam -o %s -t {{args.nthread}} %s' % (outfile, infile)
		runcmd (cmd)
	else:
		bamfile = outfile
		if informat == 'sam':
			bamfile = "{{job.outdir}}/{{in.infile | fn}}.s2b.bam"
			cmd = '{{args.sambamba}} view -S -f bam -o %s -t {{args.nthread}} %s' % (bamfile, infile)
			runcmd (cmd)
			infile = bamfile
		if doSort:
			{% if args.sortby | lambda x: x == 'queryname' %}
			params['n'] = True
			params['N'] = True
			{% endif %}
			bamfile = "{{job.outdir}}/{{in.infile | fn}}.sorted.bam"
			params['m'] = {{args.mem | quote}}
			params['tmpdir'] = tmpdir
			params['o'] = bamfile
			params['t'] = {{args.nthread}}
			cmd = '{{args.sambamba}} sort %s "%s"' % (params2CmdArgs(params, noq = ['t']), infile)
			runcmd (cmd)
			if infile != {{in.infile | quote}}:
				remove (infile)
			infile = bamfile
		if doMarkdup:
			rmdup = ""
			if doRmdup:
				rmdup = "-r"
			bamfile = "{{job.outdir}}/{{in.infile | fn}}.dedup.bam"
			cmd = '{{args.sambamba}} markdup %s -t {{args.nthread}} --tmpdir="%s" "%s" "%s"' % (rmdup, tmpdir, infile, bamfile)
			runcmd (cmd)
			if infile != {{in.infile | quote}}:
				remove (infile)
			infile = bamfile
		if doIndex:
			if path.exists (infile + '.bai'):
				move (infile + '.bai', {{out.idxfile | quote}})
			else:
				cmd = '{{args.sambamba}} index -t {{args.nthread}} "%s" "%s"' % (infile, {{out.idxfile | quote}})
				runcmd (cmd)
		if infile != outfile:
			if path.exists(infile + '.bai'):
				move (infile + '.bai', outfile + '.bai')
			move (infile, outfile)

	############### samtools
	{% elif args.tool | lambda x: x == 'samtools' %}
	if not (doSort or doIndex or doMarkdup or doRmdup):
		cmd = '{{args.samtools}} view -b -o "%s" -O bam "%s"' % (outfile, infile)
		runcmd (cmd)
	else:
		bamfile = outfile
		if doSort:
			mem = mem2({{ args.mem | quote }}, 'M')
			sortby = ''
			if {{args.sortby | quote}} == 'queryname':
				sortby = '-n'
			bamfile = "{{job.outdir}}/{{in.infile | fn}}.sorted.bam" 
			cmd = '{{args.samtools}} sort -m %sM %s -o "%s" -T "%s" -@ {{args.nthread}} -O bam "%s"' % (mem, sortby, bamfile, tmpdir, infile)
			runcmd (cmd)
			if infile != {{in.infile | quote}}:
				remove (infile)
			infile = bamfile
		if doMarkdup or doRmdup:
			bamfile = "{{job.outdir}}/{{in.infile | fn}}.dedup.bam"
			cmd = '{{args.samtools}} rmdup "%s" "%s"' % (infile, bamfile)
			runcmd (cmd)
			if infile != {{in.infile | quote}}:
				remove (infile)
			infile = bamfile
		if doIndex:
			cmd = '{{args.samtools}} index "%s" "%s"' % (bamfile, {{out.idxfile | quote}})
			runcmd (cmd)
		if infile != outfile:
			if path.exists(infile + '.bai'):
				move (infile + '.bai', outfile + '.bai')
			move (infile, outfile)

	############### picard
	{% elif args.tool | lambda x: x == 'picard' %}
	mem = mem2({{ args.mem | quote }}, 'java')
	if not (doSort or doIndex or doMarkdup or doRmdup):
		cmd = '{{args.picard}} SamFormatConverter %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="%s" O="%s"' % (mem, tmpdir, tmpdir, infile, outfile)
		runcmd (cmd)
	else:
		bamfile = outfile
		if doSort:
			bamfile = "{{job.outdir}}/{{in.infile | fn}}.sorted.bam" 
			cmd = '{{args.picard}} SortSam %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="%s" O="%s" SO={{args.sortby}}' % (mem, tmpdir, tmpdir, infile, bamfile)
			runcmd (cmd)
			if infile != {{in.infile | quote}}:
				remove (infile)
			infile = bamfile
		if doMarkdup:
			rmdup = ""
			if doRmdup:
				rmdup = "REMOVE_DUPLICATES=true"
			mfile = "/dev/null"
			bamfile = "{{job.outdir}}/{{in.infile | fn}}.dedup.bam"
			cmd = '{{args.picard}} MarkDuplicates %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="%s" O="%s" M="%s" ' % (mem, tmpdir, tmpdir, infile, bamfile, mfile)
			runcmd (cmd)
			if infile != {{in.infile | quote}}:
				remove (infile)
			infile = bamfile
		if doIndex:
			cmd = '{{args.picard}} BuildBamIndex %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="%s" O="%s"' % (mem, tmpdir, tmpdir, infile, {{out.idxfile | quote}})
			runcmd (cmd)
		if infile != outfile:
			if path.exists(infile + '.bai'):
				move (infile + '.bai', outfile + '.bai')
			move (infile, outfile)
	{% endif %}

	if not path.exists ({{out.idxfile | quote}}):
		symlink (outfile, {{out.idxfile | quote}})
		
except Exception as ex:		
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)