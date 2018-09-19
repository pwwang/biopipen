from shutil import move, rmtree
from os import makedirs, path, symlink, remove
from sys import stderr
from pyppl import Box
from bioprocs.utils import runcmd, mem2, cmdargs

infile    = {{ i.infile | quote }}
outfile   = {{ o.outfile | quote }}
informat  = {{ args.infmt | quote }}
informat  = informat if informat else {{ i.infile | ext | [1:] | quote}}
tmpdir    = path.join ("{{args.tmpdir}}", "{{proc.id}}.{{i.infile | fn}}.{{job.index}}")
doSort    = {{ args.sort }}
doIndex   = {{ args.index }}
doMarkdup = {{ args.markdup }}
doRmdup   = {{ args.rmdup }}
params    = {{ args.params }}
if doRmdup:
	doMarkdup = True
if not path.exists (tmpdir):
	makedirs (tmpdir)
try:
{% case args.tool %}
	############### biobambam
	{% when 'biobambam' %}
	mem = mem2({{ args.mem | quote }}, 'M')
	if doIndex:
		params['index'] = 1
		params['indexfilename'] = {{o.idxfile | quote}}
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


	cmd = '{{args.biobambam_bamsort}} %s' % cmdargs(params, dash = '', equal = '=')
	runcmd (cmd)

	############### sambamba
	{% when 'sambamba' %}
	if not (doSort or doIndex or doMarkdup or doRmdup):
		cmd = '{{args.sambamba}} view -S -f bam -o %s -t {{args.nthread}} %s' % (outfile, infile)
		runcmd (cmd)
	else:
		bamfile = outfile
		if informat == 'sam':
			bamfile = "{{job.outdir}}/{{i.infile | fn}}.s2b.bam"
			cmd = '{{args.sambamba}} view -S -f bam -o %s -t {{args.nthread}} %s' % (bamfile, infile)
			runcmd (cmd)
			infile = bamfile
		if doSort:
			{% if args.sortby == 'queryname' %}
			params['n'] = True
			params['N'] = True
			{% endif %}
			bamfile = "{{job.outdir}}/{{i.infile | fn}}.sorted.bam"
			params['m'] = {{args.mem | quote}}
			params['tmpdir'] = tmpdir
			params['o'] = bamfile
			params['t'] = {{args.nthread}}
			cmd = '{{args.sambamba}} sort %s "%s"' % (cmdargs(params), infile)
			runcmd (cmd)
			if infile != {{i.infile | quote}}:
				remove (infile)
			infile = bamfile
		if doMarkdup:
			rmdup = ""
			if doRmdup:
				rmdup = "-r"
			bamfile = "{{job.outdir}}/{{i.infile | fn}}.dedup.bam"
			cmd = '{{args.sambamba}} markdup %s -t {{args.nthread}} --tmpdir="%s" "%s" "%s"' % (rmdup, tmpdir, infile, bamfile)
			runcmd (cmd)
			if infile != {{i.infile | quote}}:
				remove (infile)
			infile = bamfile
		if doIndex:
			if path.exists (infile + '.bai'):
				move (infile + '.bai', {{o.idxfile | quote}})
			else:
				cmd = '{{args.sambamba}} index -t {{args.nthread}} "%s" "%s"' % (infile, {{o.idxfile | quote}})
				runcmd (cmd)
		if infile != outfile:
			if path.exists(infile + '.bai'):
				move (infile + '.bai', outfile + '.bai')
			move (infile, outfile)

	############### samtools
	{% when 'samtools' %}
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
			bamfile = "{{job.outdir}}/{{i.infile | fn}}.sorted.bam"
			cmd = '{{args.samtools}} sort -m %sM %s -o "%s" -T "%s" -@ {{args.nthread}} -O bam "%s"' % (mem, sortby, bamfile, tmpdir, infile)
			runcmd (cmd)
			if infile != {{i.infile | quote}}:
				remove (infile)
			infile = bamfile
		if doMarkdup or doRmdup:
			bamfile = "{{job.outdir}}/{{i.infile | fn}}.dedup.bam"
			cmd = '{{args.samtools}} rmdup "%s" "%s"' % (infile, bamfile)
			runcmd (cmd)
			if infile != {{i.infile | quote}}:
				remove (infile)
			infile = bamfile
		if doIndex:
			cmd = '{{args.samtools}} index "%s" "%s"' % (bamfile, {{o.idxfile | quote}})
			runcmd (cmd)
		if infile != outfile:
			if path.exists(infile + '.bai'):
				move (infile + '.bai', outfile + '.bai')
			move (infile, outfile)

	############### picard
	{% when 'picard' %}
	mem = mem2({{ args.mem | quote }}, 'java')
	if not (doSort or doIndex or doMarkdup or doRmdup):
		cmd = '{{args.picard}} SamFormatConverter %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="%s" O="%s"' % (mem, tmpdir, tmpdir, infile, outfile)
		runcmd (cmd)
	else:
		bamfile = outfile
		if doSort:
			bamfile = "{{job.outdir}}/{{i.infile | fn}}.sorted.bam"
			cmd = '{{args.picard}} SortSam %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="%s" O="%s" SO={{args.sortby}}' % (mem, tmpdir, tmpdir, infile, bamfile)
			runcmd (cmd)
			if infile != {{i.infile | quote}}:
				remove (infile)
			infile = bamfile
		if doMarkdup:
			rmdup = ""
			if doRmdup:
				rmdup = "REMOVE_DUPLICATES=true"
			mfile = "/dev/null"
			bamfile = "{{job.outdir}}/{{i.infile | fn}}.dedup.bam"
			cmd = '{{args.picard}} MarkDuplicates %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="%s" O="%s" M="%s" ' % (mem, tmpdir, tmpdir, infile, bamfile, mfile)
			runcmd (cmd)
			if infile != {{i.infile | quote}}:
				remove (infile)
			infile = bamfile
		if doIndex:
			cmd = '{{args.picard}} BuildBamIndex %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="%s" O="%s"' % (mem, tmpdir, tmpdir, infile, {{o.idxfile | quote}})
			runcmd (cmd)
		if infile != outfile:
			if path.exists(infile + '.bai'):
				move (infile + '.bai', outfile + '.bai')
			move (infile, outfile)
{% endcase %}

	if not path.exists ({{o.idxfile | quote}}):
		symlink (outfile, {{o.idxfile | quote}})

except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
